#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Luis Bonah
# Description: Loomis-Wood Plot Software for Assigning experimental Spectra to Quantum Numbers

CREDITSSTRING = """Made by Luis Bonah

As this programs GUI is based on PyQt6, which is GNU GPL v3 licensed, this program is also licensed under GNU GPL v3 (See the bottom paragraph).

pandas, matplotlib, scipy and numpy were used for this program, speeding up the development process massively.

Copyright (C) 2020

	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

##
## Global Constants and Imports
##
APP_TAG = "LLWP"

import os
import sys
import re
import io
import time
import wrapt
import random
import json
import queue
import threading
import configparser
import traceback as tb
import numpy as np
import pandas as pd
import subprocess
import webbrowser
import pyckett

from scipy import optimize, special, signal

from PyQt6.QtCore import *
from PyQt6.QtWidgets import *
from PyQt6.QtGui import *

import matplotlib
from matplotlib import style, figure
from matplotlib.backends.backend_qtagg import FigureCanvas, NavigationToolbar2QT

import warnings
try:
	warnings.simplefilter('ignore', np.RankWarning)
except AttributeError:
	warnings.simplefilter('ignore', np.exceptions.RankWarning)

QLocale.setDefault(QLocale("en_EN"))

##
## Global Decorators
##
def stopwatch_d(func):
	def timed(*args, **kwargs):
		start_time = time.time()
		result = func(*args, **kwargs)
		stop_time = time.time()
		print(f"Executing {func.__name__} took {stop_time-start_time}.")
		return result
	return timed

def askfilesfirst_d(func):
	def tmp(*args, **kwargs):
		if kwargs.get("add_files") == True:
			kwargs["add_files"] = QFileDialog.getOpenFileNames(None, f'Choose {args[1].capitalize()} File(s)',)[0]
			if len(kwargs["add_files"]) == 0:
				return
		return(func(*args, **kwargs))
	return tmp

def threading_d(func):
	def run(*args, **kwargs):
		t = threading.Thread(target=func, args=args, kwargs=kwargs)
		t.start()
		return t
	return run

def synchronized_d(lock):
	@wrapt.decorator
	def _wrapper(wrapped, instance, args, kwargs):
		with lock:
			return wrapped(*args, **kwargs)
	return _wrapper

def working_d(func):
	def wrapper(self, *args, **kwargs):
		queue_ = main.plotwidget.working
		queue_.put(1)
		if not queue_.empty():
			main.signalclass.setindicator.emit("<span style='font-weight: 600;'>Working...</span>")

		try:
			return(func(self, *args, **kwargs))
		except Exception as E:
			raise
		finally:
			queue_.get()
			queue_.task_done()
			if queue_.empty():
				main.signalclass.setindicator.emit("Ready")
	return(wrapper)

locks = {key: threading.RLock() for key in (
  "exp_df", "cat_df", "lin_df", "ser_df", "windows", "axs", "pipe", "currThread", "fitting", "peaks", "calibration",
)}




##
## Main Class, Window and Widget
##
class Main():
	def __init__(self):
		self.app = QApplication(sys.argv)
		self.app.setStyle("Fusion")
	
		self.open_windows = {}
		self.signalclass = SignalClass()
		self.config = Config(self.signalclass.updateconfig)

		self.setup_dfs()
		self.loadoptions()

	def gui(self):
		sys.excepthook = except_hook
		threading.excepthook = lambda args: except_hook(*args[:3])

		self.update_mpl_theme()
		self.app.styleHints().colorSchemeChanged.connect(self.update_mpl_theme)
		matplotlib.rcParams['axes.facecolor'] = '#00000000'
		
		matplotlib.rc('font', **self.config["plot_font_dict"])
		self.config.register("plot_font_dict", lambda: matplotlib.rc('font', **self.config["plot_font_dict"]))

		self.notificationsbox = NotificationsBox()
		self.signalclass.notification.connect(lambda text: self.notificationsbox.add_message(text))

		self.mainwindow = MainWindow()
		self.plotwidget = PlotWidget()
		self.plotwidget.gui()
		self.plotwidget.create_plots()
		
		self.mainwindow.setCentralWidget(self.plotwidget)
		self.mainwindow.loaddockables()
		self.mainwindow.createmenu()
		self.mainwindow.readstate()
		if len(sys.argv) > 1:
			project = sys.argv[1]
			if project == "RESTART":
				self.loadproject(llwpfile(".llwp"))
			else:
				self.loadproject(project)
		self.mainwindow.show()

		if self.messages:
			main.notification("\n".join(self.messages))

		self.autosavetimer = QTimer(self.app)
		self.autosavetimer.timeout.connect(lambda: main.save_lines_lin(llwpfile(".lin"), force_noappend=True, force_lin=True, quiet=True))
		if main.config["flag_autosave"]:
			self.autosavetimer.start(main.config["flag_autosave"]*1000)
		main.config.register("flag_autosave", lambda: self.autosavetimer.start(main.config["flag_autosave"]*1000) if main.config["flag_autosave"] > 0 else self.autosavetimer.stop())

		# @Luis: Add here any code that should be performed automatically after the start

		sys.exit(self.app.exec())

	def update_mpl_theme(self):
		if is_dark_theme():
			mpl_style = 'dark_background'
			mpl_background = 'black'
		else:
			mpl_style = 'default'
			mpl_background = 'white'
		
		matplotlib.style.use(mpl_style)
		matplotlib.rcParams['figure.facecolor'] = mpl_background

	def setup_dfs(self):
		self.exp_df = pd.DataFrame(columns=exp_dtypes.keys()).astype(exp_dtypes)
		self.exp_df["filename"] = None
		self.cat_df = pd.DataFrame(columns=cat_dtypes.keys()).astype(cat_dtypes)
		self.cat_df["filename"] = None
		self.lin_df = pd.DataFrame(columns=lin_dtypes.keys()).astype(lin_dtypes)
		self.lin_df["filename"] = None
		self.new_df = pd.DataFrame(columns=lin_dtypes.keys()).astype(lin_dtypes)
		self.new_df["filename"] = None

		self.ser_dtypes = {key: lin_dtypes[key] for key in list(lin_dtypes.keys())[0:13]}
		self.ser_df = pd.DataFrame(columns=self.ser_dtypes.keys()).astype(self.ser_dtypes)

		self.yrange_exp = [-1, 1]

	def notification(self, text):
		time_str = time.strftime("%H:%M", time.localtime())
		output = f"{time_str}: {text}"

		if self.config["flag_debug"] == True:
			print(output)
		if self.config["flag_shownotification"]:
			self.signalclass.notification.emit(output)
		self.signalclass.writelog.emit(output)

	def update_plot(self, dict_={}):
		self.signalclass.updateplot.emit(dict_)

	def loadoptions(self, fname=None):
		if not fname:
			self.config.update({key: value for key, (value, class_) in config_specs.items()})
			fname = llwpfile(".ini")

		config_parser = configparser.ConfigParser(interpolation=None)
		config_parser.read(fname)

		self.messages = []
		for section in config_parser.sections():
			for key, value in config_parser.items(section):
				fullkey = f"{section.lower()}_{key.lower()}"
				if fullkey in config_specs:
					try:
						class_ = config_specs[fullkey][1]
						if class_ in (dict, list, tuple):
							value = json.loads(value)
						elif class_ == bool:
							value = True if value in ["True", "1"] else False
						elif class_ == str:
							value = value.encode("utf-8").decode("unicode_escape")
						value = class_(value)
						self.config[fullkey] = value
					except Exception as E:
						message = f"The value for the option {fullkey} from the option file was not understood."
						self.messages.append(message)
						print(message)
				else:
					self.config[fullkey] = value
		
		# Special case changing colors for better contrast
		for key, value in self.config.items():
			if key in config_specs and config_specs[key][1] == Color:
				if is_dark_theme():
					if value in ("#000000", ):
						self.config[key] = "#ffffff"
						self.messages.append(f"Changed the color of '{key}' from black to white as it is otherwise invisible.")
				else:
					if value in ("#ffffff", ):
						self.config[key] = "#000000"
						self.messages.append(f"Changed the color of '{key}' from white to black as it is otherwise invisible.")
				
	def saveoptions(self):
		self.mainwindow.savestate()
		output_dict = {}
		for key, value in self.config.items():
			category, name = key.split("_", 1)
			category = category.capitalize()
			if category not in output_dict:
				output_dict[category] = {}
			if type(value) in (dict, list, tuple):
				value = json.dumps(value)
			elif type(value) == str:
				value = value.encode("unicode_escape").decode("utf-8")

			output_dict[category][name] = value

		del output_dict["Files"]

		config_parser = configparser.ConfigParser(interpolation=None)
		for section in output_dict:
			config_parser.add_section(section)
			for key in output_dict[section]:
				config_parser.set(section, key, str(output_dict[section][key]))

		with open(llwpfile(".ini"), "w+", encoding="utf-8") as file:
			config_parser.write(file)
		self.notification("Options were saved successfully!")

	def saveproject(self, fname=None):
		if fname == None:
			fname = QFileDialog.getSaveFileName(None, 'Choose file to save project to', "", "Project Files (*.llwp);;All Files (*)")[0]
		if not fname:
			return

		output_dict = {}
		for key in ["files_exp", "files_cat", "files_lin", "pipe_current", "pipe_commands"]:
			value = self.config[key]
			category, name = key.split("_", 1)
			category = category.capitalize()
			if category not in output_dict:
				output_dict[category] = {}
			if type(value) in (dict, list, tuple):
				output_dict[category][name] = json.dumps(value)
			else:
				output_dict[category][name] = value

		config_parser = configparser.ConfigParser(interpolation=None)
		for section in output_dict:
			config_parser.add_section(section)
			for key in output_dict[section]:
				config_parser.set(section, key, str(output_dict[section][key]))

		with open(fname, "w+", encoding="utf-8") as file:
			config_parser.write(file)
		self.notification("Project was saved successfully!")

	def loadproject(self, fname=None):
		if fname == None:
			fname = QFileDialog.getOpenFileName(None, 'Choose Project to load',"","Project Files (*.llwp);;All Files (*)")[0]
		if not fname:
			return

		self.loadoptions(fname)
		if self.messages:
			main.notification("\n".join(self.messages))

		return self.reread_files(do_QNs=True)

	@askfilesfirst_d
	@threading_d
	@working_d
	def load_file(self, type, keep_old=False, add_files=False, reread=False, skip_update=False, do_QNs=True):
		if reread == True:
			keep_old = False
			fnames = self.config[f"files_{type}"].keys()
		elif add_files != False:
			fnames = add_files
		else:
			fnames = []

		lock = locks[f"{type}_df"]

		if keep_old == False:
			with lock:
				df = self.return_df(type)
				df.drop(df.index, inplace=True)
				if reread == False:
					self.config[f"files_{type}"].clear()

		results = queue.Queue()
		config_updates = queue.Queue()
		errors = queue.Queue()

		if self.config["flag_loadfilesthreaded"]:
			threads = []
			for fname in fnames:
				t = threading.Thread(target=self.load_file_core, args=(fname, type, config_updates, results, errors, do_QNs))
				t.start()
				threads.append(t)

			for thread in threads:
				thread.join()

		else:
			for fname in fnames:
				# Try except block as otherwise files after an exception are not loaded
				try:
					self.load_file_core(fname, type, config_updates, results, errors, do_QNs)
				except Exception as E:
					pass

		with lock:
			for tmp_dict in list(config_updates.queue):
				self.config[f"files_{type}"].update(tmp_dict)

			df = self.return_df(type)
			if len(fnames) != 0:
				df = df[~df.filename.isin(fnames)]
			results.put(df)
			if type == "exp":
				self.exp_df = pd.concat(list(results.queue), ignore_index=True).sort_values("x", kind="merge")
			elif type == "cat":
				self.cat_df = pd.concat(list(results.queue), ignore_index=True).sort_values("x", kind="merge")
			elif type == "lin":
				self.lin_df = pd.concat(list(results.queue), ignore_index=True).sort_values("x", kind="merge")

		df = self.return_df(type)
		if type == "exp":
			self.yrange_exp = np.array((df["y"].min(), df["y"].max()))

		elif type == "cat":
			if do_QNs and len(df) and self.config["flag_autosetqns"]:
				qn_labels = ['qnu1', 'qnu2', 'qnu3', 'qnu4', 'qnu5', 'qnu6']
				QNs = len(qn_labels)
				for i in range(len(qn_labels)):
					tmp_label = qn_labels[i]
					tmp_unique = df[tmp_label].unique()
					if len(tmp_unique) == 1 and tmp_unique[0] == pyckett.SENTINEL:
						QNs = i
						break

				self.config["series_qns"] = QNs
				self.notification(f"After analysing your cat files the number of QNs was set to {QNs}.")

		elif type == "lin":
			pass

		errors = list(errors.queue)
		self.signalclass.fileschanged.emit()
		if skip_update != True:
			self.plotwidget.set_data()
		if len(fnames):
			error_text = f"<span style='color:#ff0000;'>ERROR</span>: Reading {type.capitalize()} files not successful. " if len(errors) != 0 else ''
			self.notification(f"{error_text}Read {str(len(fnames)-len(errors))+'/' if len(errors) != 0 else ''}{len(fnames)} {type.capitalize()} files successfully.")

	def load_file_core(self, fname, type, config_updates, results, errors, do_QNs):
		try:
			if not os.path.isfile(fname):
				errors.put(fname)
				self.notification(f"<span style='color:#ff0000;'>ERROR</span>: The file {fname} could not be found. Please check the file.")
			if os.path.getsize(fname) == 0:
				self.notification(f"<span style='color:#eda711;'>WARNING</span>: The file {fname} is empty and was therefore skipped.")
				return

			options = self.config[f"files_{type}"].get(fname, {})
			extension = os.path.splitext(fname)[1]

			if options.get("color") == None:
				options["color"] = self.config[f"color_{type}"]

			if type == "exp":
				args = (chr(self.config["flag_separator"]), self.config["flag_xcolumn"], self.config["flag_ycolumn"], False)
				data = exp_to_df(fname, *args)
				options["xrange"] = [data["x"].min(), data["x"].max()]
			elif type == "cat":
				formats = self.config["flag_predictionformats"]
				if extension in formats.keys():
					format = formats[extension].copy()
					intens_log = format.get("intensity_log")
					if intens_log:
						del format["intensity_log"]
					data = pd.read_fwf(fname, **format)
					data["filename"] = fname
					if intens_log:
						data["y"] = 10 ** data["y"]
					for column in cat_dtypes.keys():
						if column not in data.columns:
							data[column] = pyckett.SENTINEL
					data = data[cat_dtypes.keys()]
				else:
					data = pyckett.cat_to_df(fname, False)
			elif type == "lin":
				formats = self.config["flag_assignmentformats"]
				if extension in formats.keys():
					format = formats[extension].copy()
					data = pd.read_fwf(fname, dtype=dtypes_dict, **format)
					for column in lin_dtypes.keys():
						if column not in data.columns:
							data[column] = pyckett.SENTINEL
					data = data[lin_dtypes.keys()]
				else:
					data = pyckett.lin_to_df(fname, False)

			config_updates.put({fname: options})
			results.put(data)
		except Exception as E:
			self.notification(f"<span style='color:#ff0000;'>ERROR</span>: There occurred an error when loading the {type.capitalize()} File {fname}. Please check the file.")
			if self.config["flag_debug"]:
				tb.print_exc()
			errors.put(fname)
			raise

	@threading_d
	def reread_files(self, do_QNs=False):
		kwargs = {"reread": True, "skip_update": True, "do_QNs": do_QNs}
		threads = []
		for type in ("exp", "cat", "lin"):
			threads.append(self.load_file(type, **kwargs))

		for thread in threads:
			thread.join()
		self.plotwidget.set_data()

	@synchronized_d(locks["lin_df"])
	def save_lines_lin(self, path = None, force_noappend=False, force_lin=False, quiet=False):
		append = self.config["flag_appendonsave"]
		options = {"options": QFileDialog.Option.DontConfirmOverwrite} if append else {}

		if not path:
			path, ext = QFileDialog.getSaveFileName(None, 'Save file', '', **options)
			if not path:
				return

		catalog = self.new_df
		output = []

		handle = "a+" if append and  not force_noappend else "w+"

		with open(path, handle, encoding="utf-8") as file:
			custom_format = self.config["flag_assignmentsavefmt"]
			if custom_format and not force_lin:
				np.savetxt(file, catalog[custom_format["names"]], delimiter=custom_format.get("delimiter", " "), fmt=custom_format.get("format", '%.18e'))
			else:
				file.write(pyckett.df_to_lin(catalog))
		if not quiet:
			self.notification(f"The {len(catalog)} newly assigned lines were saved to the file {path}.")

	@synchronized_d(locks["pipe"])
	def run_pipe(self, index=None):
		if index == None:
			index = self.config["pipe_current"]

		if len(self.config["pipe_commands"]) == 0:
			self.notification("<span style='color:#eda711;'>WARNING</span>: No Pipe command specified, therefore no Pipe process was started.")
			return

		title, command, exp_rr, cat_rr, lin_rr = self.config["pipe_commands"][index]

		if command == None:
			self.notification("<span style='color:#eda711;'>WARNING</span>: No Pipe command specified, therefore no Pipe process was started.")
			return

		command = command.replace("\n", " && ")
		try:
			output = subprocess.check_output(command, shell=True)
			output = output.decode("utf-8")

			self.notification(f"The subprocess was started and returned:\n{output}")

			for type, reread in zip(("exp", "cat", "lin"), (exp_rr, cat_rr, lin_rr)):
				if reread:
					self.load_file(type, reread=True, do_QNs=False)

		except Exception as E:
			self.notification(f"The command '{command}' failed with the Exception '{E}'.")
			raise

	def get_visible_data(self, type, xrange=None, binning=False, force_all=False, scale=True):
		if type == "exp":
			with locks["exp_df"]:
				dataframe = self.exp_df.copy()
				fd = self.config["files_exp"]
		elif type == "cat":
			with locks["cat_df"]:
				dataframe = self.cat_df.copy()
				fd = self.config["files_cat"]
		elif type == "lin":
			with locks["lin_df"]:
				self.new_df["filename"] = "__lin__"
				dataframe = pd.concat((self.lin_df, self.new_df), join="inner", ignore_index=True).sort_values("x")
				fd = self.config["files_lin"]

		if xrange != None:
			x_start = dataframe["x"].searchsorted(xrange[0], side="left")
			x_stop  = dataframe["x"].searchsorted(xrange[1], side="right")
			dataframe = dataframe.iloc[x_start:x_stop].copy()

		if force_all != True:
			visible_files = {file for file in fd.keys() if not fd[file].get("hidden", False)}
			# Special Case Hide/Show catalog files
			if type == "lin" and self.config["flag_hidecatalog"] == False:
				visible_files.add("__lin__")

			if len(visible_files) != len(fd) + (type == "lin"):
				# Keep the inplace, as otherwise SettingWithCopyWarning is raised
				dataframe.query("filename in @visible_files", inplace=True)

		if binning:
			bins = main.config["plot_bins"]
			nobinning = main.config["plot_skipbinning"]
			binwidth = (xrange[1]-xrange[0]) / bins

			if len(dataframe) > max(bins, nobinning)  and binwidth != 0:
				dataframe = bin_data(dataframe, binwidth, xrange)

		if scale and type in ["exp", "cat"]:
			scalingfactordict = {file: self.config[f"files_{type}"][file].get("scale", 1) for file in fd.keys()}
			dataframe["y"] *= dataframe["filename"].replace(scalingfactordict)

		return(dataframe)

	def return_df(self, type):
		with locks[f"{type}_df"]:
			if type == "exp":
				return(self.exp_df)
			elif type == "cat":
				return(self.cat_df)
			elif type == "lin":
				return(self.lin_df)
			elif type == "new":
				return(self.new_df)
			elif type == "ser":
				return(self.ser_df)


	@synchronized_d(locks["windows"])
	def open_window(self, window, *args):
		if window not in self.open_windows:
			self.open_windows[window] = window(window.__name__.lower(), *args)
		self.open_windows[window].show()
		self.open_windows[window].activateWindow()
		return(self.open_windows[window])

	@synchronized_d(locks["lin_df"])
	def catalog_table_delete(self, all=False):
		table = main.mainwindow.catalogwindow.catalogTable
		model = main.mainwindow.catalogwindow.catalogModel
		if all:
			self.new_df.drop(self.new_df.index, inplace=True)
		else:
			selected = [index.row() for index in table.selectionModel().selectedRows()]
			self.new_df.drop(selected, inplace=True)
			self.new_df.reset_index(inplace=True, drop=True)

		self.plotwidget.set_data()
		model.update()

class MainWindow(QMainWindow):
	def __init__(self, parent=None):
		super().__init__(parent)
		self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
		self.setWindowTitle(APP_TAG)
		self.setAcceptDrops(True)
		self.shortcuts()

		geometry = main.config.get("windowgeometry_mainwindow")
		if geometry:
			self.setGeometry(*json.loads(geometry))
			
			screen_box = self.screen().geometry()
			widget_top_left = self.geometry().topLeft()
			widget_bottom_right = self.geometry().bottomRight()
			
			if not (screen_box.contains(widget_top_left) and screen_box.contains(widget_bottom_right)):
				primary_screen = QApplication.instance().primaryScreen()
				self.move(primary_screen.geometry().center()- self.rect().center())


		try:
			possible_folders = [os.path.dirname(os.path.realpath(__file__)), os.getcwd()]
			for folder in possible_folders:
				iconpath = os.path.join(folder, "LLWP.svg")
				if os.path.isfile(iconpath):
					icon = QIcon(iconpath)
					break
				
			main.app.setWindowIcon(QIcon(iconpath))
			import ctypes
			ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(APP_TAG)
		except Exception as E:
			pass

	def togglefullscreen(self):
		if self.isFullScreen():
			self.showNormal()
		else:
			self.showFullScreen()

	def savestate(self):
		state = self.saveState().data()
		geometry = self.saveGeometry().data()
		main.config["mainwindow_state"] = encodebytes(state)
		main.config["mainwindow_geometry"] = encodebytes(geometry)

	def readstate(self):
		try:
			state = main.config.get("mainwindow_state")
			geometry = main.config.get("mainwindow_geometry")
			if geometry:
				self.restoreGeometry(decodebytes(geometry))
			if state:
				self.restoreState(decodebytes(state))
		except Exception as E:
			main.notification("Could not load the saved state of the mainwindow.")

	def loaddockables(self):
		self.configureplotswindow = ConfigurePlotsWindow()
		self.referenceserieswindow = ReferenceSeriesWindow()
		self.catalogwindow = CatalogWindow()
		self.logwindow = LogWindow()
		self.quotewindow = QuoteWindow()
		self.hoverwindow = HoverWindow()
		self.hoverwindow.setVisible(False)

	def createmenu(self):
		menus = {label: self.menuBar().addMenu(f"&{label}") for label in ("Files", "View", "Fit", "Plot", "Modules", "Info")}
		for menu in menus.values():
			menu.setToolTipsVisible(True)


		##
		## Special Case Fit Function Menu
		fitfunction_menu = QMenu("Choose Fit Function", parent=self)
		self.fitfunction_actions = {}

		for function in main.plotwidget.fitfunctions:
			self.fitfunction_actions[function] = QQ(QAction, parent=self, text=f"&{function}", change=lambda x, function=function: main.config.__setitem__("fit_function", function), checkable=True, value=(main.config["fit_function"] == function))
			fitfunction_menu.addAction(self.fitfunction_actions[function])
		main.config.register("fit_function", lambda: self.change_fitfunction())

		toggleaction_configureplots = self.configureplotswindow.toggleViewAction()
		toggleaction_configureplots.setShortcut("Shift+2")
		toggleaction_configureplots.setToolTip("Toggle the visibility of the Configure Plots window")
		toggleaction_referenceseries = self.referenceserieswindow.toggleViewAction()
		toggleaction_referenceseries.setShortcut("Shift+3")
		toggleaction_referenceseries.setToolTip("Toggle the visibility of the Reference Series window")
		toggleaction_catalog = self.catalogwindow.toggleViewAction()
		toggleaction_catalog.setShortcut("Shift+4")
		toggleaction_catalog.setToolTip("Toggle the visibility of the Catalog of newly assigned Lines window")
		toggleaction_log = self.logwindow.toggleViewAction()
		toggleaction_log.setShortcut("Shift+5")
		toggleaction_log.setToolTip("Toggle the visibility of the Log window")
		toggleaction_hover = self.hoverwindow.toggleViewAction()
		toggleaction_hover.setShortcut("Shift+6")
		toggleaction_hover.setToolTip("Toggle the visibility of the Hover window")

		actions_to_menus = {
			"Files": (
				QQ(QAction, parent=self, text="&Load Spectrum", change=lambda x: main.load_file("exp", add_files=True), tooltip="Replace Exp file(s)"),
				QQ(QAction, parent=self, text="&Add Spectrum", change=lambda x: main.load_file("exp", add_files=True, keep_old=True), tooltip="Add Exp file(s)"),
				QQ(QAction, parent=self, text="&Load Cat File", change=lambda x: main.load_file("cat", add_files=True), tooltip="Replace Cat file(s)"),
				QQ(QAction, parent=self, text="&Add Cat File", change=lambda x: main.load_file("cat", add_files=True, keep_old=True), tooltip="Add Cat file(s)"),
				QQ(QAction, parent=self, text="&Load Lin File", change=lambda x: main.load_file("lin", add_files=True), tooltip="Replace Lin file(s)"),
				QQ(QAction, parent=self, text="&Add Lin File", change=lambda x: main.load_file("lin", add_files=True, keep_old=True), tooltip="Add Lin file(s)"),
				None,
				QQ(QAction, parent=self, text="&Reread Files", change=main.reread_files, tooltip="Reread all Exp, Cat and Lin files", shortcut="Ctrl+R"),
				None,
				QQ(QAction, parent=self, text="&Edit Files", shortcut="Shift+7", tooltip="See current files and their options", change=lambda x: main.open_window(FileWindow)),
				None,
				QQ(QAction, parent=self, text="&Save current values as default", shortcut="Ctrl+D", tooltip="Save current configuration as default", change=lambda x: main.saveoptions()),
				None,
				QQ(QAction, parent=self, text="&Save Project", change=lambda x: main.saveproject(), shortcut="Ctrl+S", tooltip="Save the references to the currently opened files as a project"),
				QQ(QAction, parent=self, text="&Load Project", change=lambda x: main.loadproject(), shortcut="Ctrl+O", tooltip="Load a project file"),
				None,
				QQ(QAction, parent=self, text="&Quit", change=self.close, tooltip="Close the program"),
			),
			"View": (
				QQ(QAction, "layout_mpltoolbar", parent=self, text="&MPL Toolbar", shortcut="Shift+1", tooltip="Show or hide toolbar to edit or save the plot canvas", checkable=True),
				toggleaction_configureplots,
				toggleaction_referenceseries,
				toggleaction_catalog,
				toggleaction_hover,
				toggleaction_log,
				None,
				QQ(QAction, "flag_alwaysshowlog", parent=self,  text="&Force Show Log", tooltip="Make log window visible if a new message is shown", checkable=True),
			),
			"Fit": (
				fitfunction_menu,
				QQ(QAction, parent=self, text="&Change Function", shortcut="Ctrl+F", tooltip="Cycle through the available fit-functions", change=lambda x: self.change_fitfunction(True)),
				None,
				QQ(QAction, "fit_alwaysfit", parent=self, text="&Always Fit", tooltip="Always fit on selecting a range and never zoom", checkable = True),
				QQ(QAction, parent=self, text="&Change Fit Color", tooltip="Change the color of the fitfunction", change=lambda x: main.plotwidget.change_fitcolor()),
				None,
				QQ(QAction, "fit_alwaysassign", parent=self, text="&Always Assign", tooltip="Always assign fitted line and add to assigned lines table", checkable = True),
				QQ(QAction, parent=self, text="&Manual Assign", tooltip="Manually assign and add line to the assigned lines table, useful when always assign is not selected", change=lambda x: main.plotwidget.manual_assign(), shortcut="Ctrl+Return"),
			),
			"Plot": (
				QQ(QAction, parent=self, text="&# Plots", shortcut="Ctrl+N", tooltip="Change number of plots", change=lambda x: main.plotwidget.plot_number()),
				QQ(QAction, parent=self, text="&Go to ..", tooltip="Go to specific frequency", shortcut="Ctrl+G", change=lambda x: main.plotwidget.offset_dialog()),
				QQ(QAction, parent=self, text="&Set Width", tooltip="Set specific width", shortcut="Ctrl+W", change=lambda x: main.plotwidget.width_dialog()),
				None,
				QQ(QAction, "flag_automatic_draw", parent=self, text="&Automatic Draw", tooltip="Update canvas automatically when plot is updated, might be switched off if program is unresponsive", checkable = True),
				QQ(QAction, parent=self, text="&Manual Draw", tooltip="Draw canvas manually", change=lambda x: main.plotwidget.manual_draw(), shortcut="Shift+Space"),
			),
			"Modules": (
				QQ(QAction, parent=self, text="&Config", shortcut="Ctrl+1", tooltip="Press to open the config module", change=lambda x: main.open_window(ConfigWindow)),
				None,
				QQ(QAction, parent=self, text="&Blended Lines", shortcut="Ctrl+2", tooltip="Press to open the blended lines module, which allows to fit and assign resolvable blends", change=lambda x: main.open_window(BlendedLinesWindow)),
				QQ(QAction, parent=self, text="&Seriesfinder", shortcut="Ctrl+3", tooltip="Press to open the seriesfinder module, which allows to explore the .cat lines, e.g. finding the most intense transitions", change=lambda x: main.open_window(SeriesfinderWindow)),
				QQ(QAction, parent=self, text="&Peakfinder", shortcut="Ctrl+4", tooltip="Press to open the peakfinder module, which allows to find the peaks of your experimental spectrum and check, which peaks are not assigned yet", change=lambda x: main.open_window(PeakfinderWindow)),
				QQ(QAction, parent=self, text="&Residuals", shortcut="Ctrl+5", tooltip="Press to open the Residuals module", change=lambda x: main.open_window(ResidualsWindow)),
				QQ(QAction, parent=self, text="&Energy Levels", shortcut="Ctrl+6", tooltip="Press to open the Energy Levels module", change=lambda x: main.open_window(EnergyLevelsWindow)),
				QQ(QAction, parent=self, text="&Spectra Resolver", shortcut="Ctrl+7", tooltip="Press to open the Spectra Resolver module", change=lambda x: main.open_window(SpectraResolverWindow)),
				QQ(QAction, parent=self, text="&Lineshape", shortcut="Ctrl+8", tooltip="Press to open the lineshape module, which allows to simulate different line profiles from the .cat stick spectrum", change=lambda x: main.open_window(LineshapeWindow)),
				QQ(QAction, parent=self, text="&Series Fit", shortcut="Ctrl+9", tooltip="Open Series Fit Window, which allows to try out different combinations of transitions as a series", change=lambda x: main.open_window(SeriesFitWindow)),
				QQ(QAction, parent=self, text="&Energy Level Trends", tooltip="Open Energy Levels Trend Window, which allows to find trends in energy level series, mitigating label shifts", change=lambda x: main.open_window(EnergyLevelsTrendWindow)),
				QQ(QAction, parent=self, text="&Create Report", tooltip="Open Report Window, which allows to summarise your analysis", change=lambda x: main.open_window(ReportWindow)),
				QQ(QAction, parent=self, text="&Create Figure", tooltip="Open Figure Window, which allows to create publication quality figures", change=lambda x: main.open_window(FigureWindow)),
				QQ(QAction, parent=self, text="&Calibrate Spectrum", tooltip="Open Calibrate Spectrum Window, which allows to intensity calibrate the spectrum", change=lambda x: main.open_window(CalibrateSpectrumWindow)),
				None,
				QQ(QAction, parent=self, text="&Pipe", shortcut="Ctrl+0", tooltip="Set up or execute the pipe command", change=lambda x: main.open_window(PipeWindow)),
				QQ(QAction, parent=self, text="&Run Pipe", tooltip="Run the pipe command", change=lambda x: main.run_pipe(), shortcut="Ctrl+P"),
			),
			"Info": (
				QQ(QAction, parent=self, text="&Send Mail to Author", tooltip="Send a mail to the developer", change=lambda x: send_mail_to_author()),
				QQ(QAction, parent=self, text="&Credits and License", tooltip="See the Credits and License", change=lambda x: main.open_window(CreditsWindow)),
			),
		}

		for label, menu in menus.items():
			for widget in actions_to_menus[label]:
				if widget is None:
					menu.addSeparator()
				elif isinstance(widget, QAction):
					menu.addAction(widget)
				else:
					menu.addMenu(widget)

	def change_fitfunction(self, nextfunction=False):
		if nextfunction:
			fitfunctions = main.plotwidget.fitfunctions
			newindex = fitfunctions.index(main.config["fit_function"])+1
			value = fitfunctions[newindex%len(fitfunctions)]
			main.config["fit_function"] = value
			return

		value = main.config["fit_function"]
		for fitfunction, action in main.mainwindow.fitfunction_actions.items():
			action.setChecked(fitfunction == value)
		main.notification(f"Fitting with {value}")

	def shortcuts(self):
		shortcuts_dict = {
			"w": lambda: main.plotwidget.set_width("++"),
			"s": lambda: main.plotwidget.set_width("--"),
			"a": lambda: main.plotwidget.set_offset("--"),
			"d": lambda: main.plotwidget.set_offset("++"),

			"Shift+w": lambda: main.plotwidget.set_width("+"),
			"Shift+s": lambda: main.plotwidget.set_width("-"),
			"Shift+a": lambda: main.plotwidget.set_offset("-"),
			"Shift+d": lambda: main.plotwidget.set_offset("+"),

			"Ctrl+k": lambda: commandline(),
			"Ctrl+Shift+k": lambda: commandline(False),
			"Shift+0": lambda: main.mainwindow.quotewindow.toggle(),

			"F11": lambda: main.mainwindow.togglefullscreen(),

			"Ctrl+l": lambda: [
			  main.config.__setitem__("pipe_current", (main.config["pipe_current"]+1) % (len(main.config["pipe_commands"] or 1))),
			  main.notification(f"Set current pipe to tab '{main.config['pipe_commands'][main.config['pipe_current']][0]}'."),
			],
			"Ctrl+Shift+l": lambda: [
			  main.config.__setitem__("commandlinedialog_current", (main.config["commandlinedialog_current"]+1) % len(main.config["commandlinedialog_commands"])),
			  main.notification(f"Set current commandline command to tab '{main.config['commandlinedialog_commands'][main.config['commandlinedialog_current']][0]}'."),
			],
		}


		for keys, function in shortcuts_dict.items():
			QShortcut(keys, self).activated.connect(function)

	def dragEnterEvent(self, event):
		if event.mimeData().hasUrls():
			event.accept()
		else:
			event.ignore()

	def dropEvent(self, event):
		files = [url.toLocalFile() for url in event.mimeData().urls()]
		files_dropped = {}

		types = main.config["flag_extensions"]
		files_by_class = {key: [] for key in list(types.keys())}

		for file in files:
			if not os.path.isfile(file):
				main.notification(f"<span style='color:#ff0000;'>ERROR</span>: The file {file} could not be found.")
				continue

			path, extension = os.path.splitext(file)
			extension = extension if extension else os.path.basename(path)
			type = None
			for key, value in types.items():
				if extension in value:
					type = key
					break

			if type == None:
				item, ok = QInputDialog.getItem(self, "Choose File Type", f"Choose the file type for the extension \"{extension}\":", [x.capitalize() for x in types], editable=False)
				if not (ok and item):
					continue
				types[item.lower()].append(extension)
				type = item.lower()

			files_by_class[type].append(file)
		self.dropEvent_core(files_by_class)

	@threading_d
	def dropEvent_core(self, files_by_class):
		threads = []
		for type in ["exp", "cat", "lin"]:
			files = files_by_class[type]
			if files:
				threads.append(main.load_file(type, keep_old=True, add_files=files, skip_update=True))

		for thread in threads:
			thread.join()

		for project in files_by_class["project"]:
			main.loadproject(project)

		main.plotwidget.set_data()

	def closeEvent(self, *args, **kwargs):
		if not main.new_df.empty:
			main.save_lines_lin(llwpfile(".lin"), force_noappend=True, force_lin=True)
		for window in list(main.open_windows.values()):
			window.close()
		return super().closeEvent(*args, **kwargs)

	def resizeEvent(self, event):
		main.config["windowgeometry_mainwindow"] = self.geometry().getRect()

	def moveEvent(self, event):
		main.config["windowgeometry_mainwindow"] = self.geometry().getRect()

class PlotWidget(QGroupBox):
	def __init__(self, parent=None):
		super().__init__(parent)

		self.fig = figure.Figure(dpi=main.config["plot_dpi"])
		self.cid = self.fig.canvas.mpl_connect('button_press_event', self.on_click)

		main.config.register("plot_dpi", lambda: self.fig.set_dpi(main.config["plot_dpi"]))
		self.plotcanvas = FigureCanvas(self.fig)
		self.plotcanvas.setMinimumHeight(200)
		self.plotcanvas.setMinimumWidth(200)
		self.fig.canvas.mpl_connect("motion_notify_event", self.on_hover)
		self.cache_positions = None
		main.signalclass.drawplot.connect(lambda: self.plotcanvas.draw())

		self.mpltoolbar = NavigationToolbar2QT(self.plotcanvas, self)
		self.mpltoolbar.setVisible(main.config["layout_mpltoolbar"])
		main.config.register("layout_mpltoolbar", lambda: self.mpltoolbar.setVisible(main.config["layout_mpltoolbar"]))

		toplayout = QHBoxLayout()

		buttonsdict = {
			"in":			lambda x: self.set_width("++"),
			"out":			lambda x: self.set_width("--"),
			"left":			lambda x: self.set_offset("-"),
			"right":		lambda x: self.set_offset("+"),
		}

		for label, func in buttonsdict.items():
			button = QQ(QPushButton, text=label, change=func, visible=main.config["flag_showmainplotcontrols"])
			toplayout.addWidget(button)
			main.config.register("flag_showmainplotcontrols", lambda button=button: button.setVisible(main.config["flag_showmainplotcontrols"]))

		self.toplabel = QQ(QLabel, text="", wordwrap=False)
		self.indicator = QQ(QLabel, text="Ready", textFormat=Qt.TextFormat.RichText)
		self.working = queue.Queue()
		main.signalclass.setindicator.connect(self.indicator.setText)

		toplayout.addWidget(self.toplabel, 1)
		toplayout.addWidget(self.indicator)

		rows_cols_elements = (
			QQ(QLabel, text="||  Plots: ", visible=main.config["flag_showmainplotrowscols"]),
			QQ(QSpinBox, "plot_rows", range=(1, None), maxWidth=45, visible=main.config["flag_showmainplotrowscols"]),
			QQ(QLabel, text="x", visible=main.config["flag_showmainplotrowscols"]),
			QQ(QSpinBox, "plot_cols", range=(1, None), maxWidth=45, visible=main.config["flag_showmainplotrowscols"]),
		)
		for elem in rows_cols_elements:
			toplayout.addWidget(elem)
			main.config.register("flag_showmainplotrowscols", lambda elem=elem: elem.setVisible(main.config["flag_showmainplotrowscols"]))

		width_elements = (
			QQ(QLabel, text="    Width: ", visible=main.config["flag_showmainplotwidth"]),
			QQ(QDoubleSpinBox, "plot_width", range=(0, None), minWidth=85, visible=main.config["flag_showmainplotwidth"], change=lambda: main.config["plot_coupled"] and self.set_data()),
		)
		for elem in width_elements:
			toplayout.addWidget(elem)
			main.config.register("flag_showmainplotwidth", lambda elem=elem: elem.setVisible(main.config["flag_showmainplotwidth"]))

		layout = QVBoxLayout()
		layout.addLayout(toplayout)
		layout.addWidget(self.plotcanvas, 1)
		layout.addWidget(self.mpltoolbar)
		self.setLayout(layout)

		self.lcp = (0, 0)
		self.lastassignment = None
		self.create_plots_id = None
		self.set_data_id = None
		self.fitfunctions = ("Pgopher", "Polynom", "Polynom Autorank", "Gauss", "Lorentz", "Voigt", "Gauss 1st Derivative", "Lorentz 1st Derivative", "Voigt 1st Derivative", "Gauss 2nd Derivative", "Lorentz 2nd Derivative", "Voigt 2nd Derivative")
		self.fitcurve = None
		self.fitline = None

		main.signalclass.updateplot.connect(lambda: self.set_data())
		main.config.register(("series_annotate_xs", "plot_annotate", "plot_bins", "flag_automatic_draw"), lambda: self.set_data())

	def gui(self):
		self.create_plots().join()
		main.config.register(("plot_rows", "plot_cols"), self.create_plots)
		main.config.register("series_references", self.set_data)

	def set_offset(self, value, absolute=True):
		coupled = main.config["plot_coupled"]
		lcp = self.get_current_plot()

		reference = self.get_series_reference()
		if reference["method"] == "Fortrat":
			offset = reference["fortrat"]["center"]
			width = main.config["plot_width"]
		elif coupled:
			offset = main.config["plot_offset"]
			width = main.config["plot_width"]
		else:
			offset = self.axs["offset"][lcp]
			width = self.axs["width"][lcp]

		if value == "+":
			offset += width/4
		elif value == "-":
			offset -= width/4
		elif value == "++":
			offset += width/2
		elif value == "--":
			offset -= width/2
		else:
			if reference["method"] == "Fortrat":
				offset = value
			else:
				if absolute:
					xpos  = self.cache_positions[lcp]
					value = value - xpos
				offset = value

		if reference["method"] == "Fortrat":
			reference["fortrat"]["center"] = offset
		elif coupled:
			main.config["plot_offset"] = offset
		else:
			with locks["axs"]:
				self.axs["offset"][lcp] = offset
		self.set_data()

	def get_offset(self, index=None, all=False):
		if main.config["plot_coupled"]:
			return main.config["plot_offset"]
		else:
			if all:
				return(self.axs["offset"])
			if index is None:
				index = self.get_current_plot()
			return self.axs["offset"][index]

	def reset_offsets(self):
		with locks["axs"]:
			shape = self.axs["ax"].shape
			self.axs["offset"] = np.zeros(shape)
			self.axs["width"]  = np.full(shape, main.config["plot_width"], dtype=np.float64)
		main.config["plot_offset"] = 0
		self.set_data()

	def set_width(self, value, absolute=True):
		coupled = main.config["plot_coupled"]
		if coupled:
			width = main.config["plot_width"]
		else:
			lcp = self.get_current_plot()
			width = self.axs["width"][lcp]

		if value == "+":
			width *= 3/4
		elif value == "-":
			width /= 3/4
		elif value == "++":
			width *= 1/2
		elif value == "--":
			width /= 1/2
		elif absolute:
			width = value
		else:
			width *= value

		if coupled:
			main.config["plot_width"] = width
		else:
			with locks["axs"]:
				self.axs["width"][lcp] = width
		self.set_data()

	def get_widths(self):
		references = main.config["series_references"]
		is_fortrat = [reference["method"] == "Fortrat" for reference in references]

		if main.config["plot_coupled"]:
			widths = np.full(self.axs["ax"].shape, main.config["plot_width"])
		else:
			widths = self.axs["width"]

		for i, do_fortrat in zip(range(widths.shape[1]), is_fortrat):
			if do_fortrat:
				reference = references[i]
				jstart, center = reference["fortrat"]["jstart"], reference["fortrat"]["center"]
				widths[:, i] *= (np.arange(0, len(widths[:, i]))[::-1] + jstart)*2

		return widths

	def get_positions(self, return_qns=False, cat_df=None):
		shape = self.axs["ax"].shape
		positions = np.zeros(shape)
		if return_qns:
			allqns = np.empty(shape, dtype=object)
		references = main.config["series_references"]

		for i, reference in zip(range(shape[1]), references):
			method = reference["method"]
			if method == "Transition":
				qns, qnus, qnls, diffs = self.get_qns(reference["transition"], return_all=True)
				if return_qns:
					allqns[:, i] = qns
				file = reference["transition"].get("file")
				positions[:, i] = self.get_position_from_qns(qns, qnus, qnls, diffs, file=file, cat_df=cat_df)
			elif method == "List":
				i0 = reference["list"]["i0"]
				xs_all = reference["list"]["xs"]

				xs = np.zeros(shape[0])
				if xs_all is not None and i0 < len(xs_all):
					imax = min(len(xs), len(xs_all)-i0)
					xs[:imax] = xs_all[i0:imax+i0]
				else:
					imax = 0
				positions[:, i] = xs[::-1]
				
				if return_qns:
					qns = reference["list"]["qns"]
					if qns is not None:
						tmp = [([], []) for _ in range(shape[0])]
						
						for j in range(imax):
							tmp[j] = qns[i0+j]
						allqns[:, i] = tmp[::-1]
					else:
						allqns[:, i] = None
				
			elif method == "Expression":
				if return_qns:
					allqns[:, i] = None
				N0 = reference["expression"]["N0"]
				try:
					expr = reference["expression"]["expression"]
					if expr:
						xs = [eval(expr, {"N": i, "N0": N0}) for i in range(shape[0])]
						positions[:, i] = xs[::-1]
					else:
						positions[:, i] = 0
				except Exception as E:
					positions[:, i] = 0
					main.notification("<span style='color:#eda711;'>WARNING</span>: Could not evaluate your expression.")
					raise
			elif method == "Fortrat":
				jstart, center = reference["fortrat"]["jstart"], reference["fortrat"]["center"]
				positions[:, i] = (np.arange(0, shape[0])[::-1] + jstart)*2*center

			self.cache_positions = positions
		if return_qns:
			return(positions, allqns)
		else:
			return(positions)

	def get_series_reference(self, column=None):
		if not column:
			column = self.get_current_plot()[1]
		references = main.config["series_references"]

		if column < len(references):
			return(references[column])
		else:
			return({"method": None})

	def get_qns_list(self, reference, index):
		qns = reference["list"]["qns"]
		if qns is None:
			return None
		
		nop = self.axs["ax"].shape[0]
		i0 = reference["list"]["i0"]
		tmp_index = nop - index[0] - 1 + i0
		if tmp_index < len(qns):
			qns = qns[tmp_index]
			return(qns)

	def get_qns(self, transition, return_all=False):
		nop = self.axs["ax"].shape[0]
		noq = main.config["series_qns"]

		qnus, qnls = np.array(transition["qnus"][:noq]), np.array(transition["qnls"][:noq])
		diffs = np.array(transition["qndiffs"][:noq]) if transition["increase_freely"] else np.array(transition["qnincrs"][:noq])

		# @Luis: Think about alternating signs
		# alternating_signs = np.array([1, -1, 1, 1, 1, 1][:noq])

		qns = []

		for i in range(nop):
			ind = nop-i-1
			# upper, lower = (qnus+ind*diffs) * alternating_signs ** ind, (qnls+ind*diffs) * alternating_signs ** ind
			upper, lower = qnus+ind*diffs, qnls+ind*diffs
			qns.append((upper, lower))

		if return_all:
			return(qns, qnus, qnls, diffs)
		else:
			return(qns)

	def get_position_from_qns(self, qns, qnus, qnls, diffs, file=None, cat_df=None):
		positions = []
		if cat_df is None:
			cat_df = main.get_visible_data("cat", scale=False)
		if file:
			cat_df = cat_df.query("filename == @file")

		# Prefiltering cat_df (Similar to seriesfitwindow -> only one big query)
		noq = main.config["series_qns"]

		conditions = []
		conditions_incr = []
		for i, qnu, qnl, diff in zip(range(noq), qnus, qnls, diffs):
			diff = int(diff)
			if diff:
				conditions_incr.append(f"((qnu{i+1} - {qnu})/{diff})")
				conditions_incr.append(f"((qnl{i+1} - {qnl})/{diff})")
			else:
				conditions.append(f"(qnu{i+1} == {qnu})")
				conditions.append(f"(qnl{i+1} == {qnl})")

		if len(conditions_incr):
			conditions.append(" == ".join(conditions_incr))

		conditions = " and ".join(conditions)
		cat_df = cat_df.query(conditions)

		for qnus, qnls in qns:
			cond_upper = [f"(qnu{i+1} == {qn})" for i, qn in enumerate(qnus)]
			cond_lower = [f"(qnl{i+1} == {qn})" for i, qn in enumerate(qnls)]
			condition  = " & ".join(cond_upper + cond_lower)

			vals = cat_df.query(condition)["x"].to_numpy()
			val = vals[0] if len(vals) else 0
			positions.append(val)
		return(positions)

	def plot_annotations(self, xpos, allqns):
		ann_dict = main.config["plot_annotations_dict"]

		if not main.config["plot_annotate"]:
			if not (self.axs["annotation"] == None).all():
				for i in range(self.axs["ax"].shape[0]):
					for j in range(self.axs["ax"].shape[1]):
						annotation = self.axs["annotation"][i, j]
						if annotation:
							annotation.remove()
							annotation.set_visible(False) # Keep this, becaue of matplotlib bug
							del annotation
						with locks["axs"]:
							self.axs["annotation"][i, j] = None
			return

		lin_df = main.get_visible_data("lin")
		for i in range(self.axs["ax"].shape[0]):
			for j in range(self.axs["ax"].shape[1]):
				annotation = self.axs["annotation"][i, j]
				x = xpos[i, j]
				qns = allqns[i, j]
				color = matplotlib.rcParams['text.color']

				if main.config["series_annotate_xs"] or qns is None:
					text = f"{{:{main.config['series_annotate_fmt']}}}".format(x)
				else:
					text = f"{', '.join([str(qn) if qn != pyckett.SENTINEL else '-' for qn in qns[0]])} ← {', '.join([str(qn) if qn != pyckett.SENTINEL else '-' for qn in qns[1]])}"

					query = " and ".join([f"qnu{i+1} == {qn}" for i, qn in enumerate(qns[0])] + [f"qnl{i+1} == {qn}" for i, qn in enumerate(qns[1])])
					if query and len(lin_df.query(query)):
						color = main.config["color_lin"]

				if not annotation:
					ax  = self.axs["ax"][i, j]
					with locks["axs"]:
						self.axs["annotation"][i, j] = ax.text(**ann_dict, s=text, transform=ax.transAxes)
				else:
					with locks["axs"]:
						self.axs["annotation"][i, j].set_text(text)
				self.axs["annotation"][i, j].set_color(color)

	def check_blends(self, index, dict_):
		with locks["axs"]:
			main.signalclass.drawplot.emit()
		if main.config["series_blenddialog"] and self.get_series_reference(index[1])["method"] == "Transition":
			blendwidth = main.config["series_blendwidth"]
			xrange = (dict_["xpre"]-blendwidth, dict_["xpre"]+blendwidth)
			entries = main.get_visible_data("cat", xrange=xrange, scale=False)
			if len(entries) > 1:
				if BlendWindow in main.open_windows:
					main.open_windows[BlendWindow].update_gui(entries, dict_, index)
					main.open_windows[BlendWindow].show()
					main.open_windows[BlendWindow].activateWindow()
				else:
					main.open_window(BlendWindow, entries, dict_, index)
				return(True)
			else:
				return(False)

		return(False)

	def on_hover(self, event):
		x = event.xdata
		y = event.ydata

		if not all([x, y, event.inaxes]):
			text_top = ""
			text_annotation = ""
		else:
			if main.config["flag_showmainplotposition"]:
				text_top = f"({x=:{main.config['flag_xformatfloat']}}, {y=:{main.config['flag_xformatfloat']}})"
			else:
				text_top = ""

			cutoff = main.config["plot_hover_cutoff"]
			xrange = (x-cutoff, x+cutoff)
			cat_df = main.get_visible_data("cat", xrange=xrange)
			lin_df = main.get_visible_data("lin", xrange=xrange)

			dataframes = {"cat": cat_df, "lin": lin_df}
			transitions = {}
			noq = main.config["series_qns"]

			for type, df in dataframes.items():
				if len(df):
					df["dist"] = abs(df["x"] - x)
					smallest_distance = df["dist"].min()
					df = df.query("dist == @smallest_distance")

					tmp = []
					for i, row in df.iterrows():
						qnus = [row[f"qnu{i+1}"] for i in range(noq)]
						qnls = [row[f"qnl{i+1}"] for i in range(noq)]
						tmp.append(f"{', '.join([str(qn) for qn in qnus if qn != pyckett.SENTINEL])} ← {', '.join([str(qn) for qn in qnls if qn != pyckett.SENTINEL])}")

					transitions[type] = tmp

			text_annotation = []
			if "cat" in transitions:
				text_annotation.append("Cat:\n" + "\n".join(transitions["cat"]))
			if "lin" in transitions:
				text_annotation.append("Lin:\n" + "\n".join(transitions["lin"]))

			if text_annotation:
				text_annotation = "\n\n".join(text_annotation)
			else:
				text_annotation = ""

		main.signalclass.writehover.emit(text_annotation)
		self.toplabel.setText(text_top)

	def on_click(self, event):
		ax = event.inaxes
		index = np.asarray(np.where(self.axs["ax"] == ax)).T
		if len(index):
			self.lcp = tuple(index[0])

	@synchronized_d(locks["axs"])
	def on_range(self, xmin, xmax, index):
		axrange = self.axs["ax"][index].get_xlim()
		if xmax == xmin or xmax > axrange[1] or xmin < axrange[0]:
			return

		shift = (QApplication.keyboardModifiers() == Qt.KeyboardModifier.ShiftModifier)
		if shift or main.config["fit_alwaysfit"]:
			xmiddle, xuncert = self.fit_peak(xmin, xmax, index)
			if main.config["fit_clipboard"]:
				QApplication.clipboard().setText(str(xmiddle))
			xpre = self.cache_positions[index]

			dict_ = {
				"x":		xmiddle,
				"error":	xuncert,
				"xpre":		xpre
			}

			reference = self.get_series_reference(index[1])
			if reference["method"] == "Transition":
				qns = self.get_qns(reference["transition"])[index[0]]
				for i, qnu, qnl in zip(range(6), qns[0], qns[1]):
					dict_[f"qnu{i+1}"] = qnu
					dict_[f"qnl{i+1}"] = qnl
			elif reference["method"] == "List":
				qns = self.get_qns_list(reference, index=index)
				if qns is not None:
					for i, qnu, qnl in zip(range(6), qns[0], qns[1]):
						dict_[f"qnu{i+1}"] = qnu
						dict_[f"qnl{i+1}"] = qnl

			if not main.config["fit_alwaysassign"]:
				self.lastassignment = index, dict_
			elif self.check_blends(index, dict_):
				return
			else:
				self.assign(index, dict_)
		else:
			center, width = (xmin + xmax) / 2, xmax - xmin
			reference = self.get_series_reference()
	
			if reference["method"] == "Fortrat":
				lcp = self.get_current_plot()
				jstart = reference["fortrat"]["jstart"]
				nplots = self.axs["ax"].shape[0]
				jplot = jstart + nplots - lcp[0] - 1
				
				center = center / jplot / 2
				width = width / jplot / 2
			
			self.set_offset(center, absolute=True)
			self.set_width(width, absolute=True)

	def fit_peak(self, xmin, xmax, index=None):
		data = main.get_visible_data("exp", xrange=(xmin, xmax))
		exp_xs = data["x"].to_numpy()
		exp_ys = data["y"].to_numpy()
		fit_xs = np.linspace(xmin, xmax, main.config["fit_xpoints"])

		amplitude_direction = main.config["fit_amplitudedirection"]
		fitfunction = main.config["fit_function"]
		if len(exp_xs) < 2 and fitfunction != "Pgopher":
			main.notification("<span style='color:#eda711;'>WARNING</span>: You did select less than two points of your spectrum, this fit will not work.")
			raise CustomError("Too few points selected")

		if self.fitline != None:
			self.fitline.remove()
			self.fitline = None
		if self.fitcurve != None:
			self.fitcurve.remove()
			self.fitcurve = None

		# @Luis: 
		try:
			if fitfunction == "Pgopher":
				ymin, ymax = np.min(exp_ys), np.max(exp_ys)
				
				if amplitude_direction < 0:
					cutoff = ymin + (ymax-ymin)/2
					mask = (exp_ys <= cutoff)
				else:
					cutoff = ymax - (ymax-ymin)/2
					mask = (exp_ys >= cutoff)

				fit_xs = exp_xs[mask]
				fit_ys = exp_ys[mask]

				# @Luis: Check for the correct sign of exp_ys here

				# @Luis: Check if average here should use relative fit_ys
				# to correct for offsets in intensity
				xmiddle = np.sum(fit_xs*fit_ys)/np.sum(fit_ys)
				xuncert = 0
			elif fitfunction == "Polynom":
				try:
					popt = np.polyfit(exp_xs, exp_ys, main.config["fit_polynomrank"])
				except Exception as E:
					popt = np.polyfit(exp_xs, exp_ys, main.config["fit_polynomrank"])
				polynom = np.poly1d(popt)
				fit_ys = polynom(fit_xs)

				if amplitude_direction < 0:
					xmiddle = fit_xs[np.argmin(fit_ys)]
				else:
					xmiddle = fit_xs[np.argmax(fit_ys)]
				
				xuncert = 0 # @Luis: Find real error for this part

			elif fitfunction == "Polynom Autorank":
				rms_best = np.inf
				rank_best = 0

				for rank in range(min((len(exp_ys), main.config["fit_polynommaxrank"]))):
					try:
						try:
							popt = np.polyfit(exp_xs, exp_ys, rank)
						except Exception as E:
							popt = np.polyfit(exp_xs, exp_ys, rank)
						polynom = np.poly1d(popt)
						fit_ys = polynom(exp_xs)

						rms = np.mean((fit_ys - exp_ys)**2)
						if rms < rms_best:
							rms_best = rms
							rank_best = rank
					except Exception as E:
						continue

				popt = np.polyfit(exp_xs, exp_ys, rank_best)
				polynom = np.poly1d(popt)
				fit_ys = polynom(fit_xs)

				if amplitude_direction < 0:
					xmiddle = fit_xs[np.argmin(fit_ys)]
				else:
					xmiddle = fit_xs[np.argmax(fit_ys)]
				xuncert = 0 # @Luis: Find real error for this part

			else:
				x0 = (xmin+xmax)/2
				ymin, ymax = np.min(exp_ys), np.max(exp_ys)
				ymean = np.mean(exp_ys)
				yptp = ymax-ymin
				w0 = (xmax-xmin)
				y0 = 0

				amp_min, amp_max = -3*yptp, 3*yptp
				if amplitude_direction < 0:
					amp_max = 0
					y0 = -yptp
				if amplitude_direction > 0:
					amp_min = 0
					y0 = yptp

				if main.config["fit_offset"]:
					function, p0, bounds = {
						"Gauss":					(lambda *x: lineshape("Gauss", 0, *x[:-1])+x[-1],   (x0, y0, w0, ymean),      ((xmin, amp_min, 0, ymin),    (xmax, amp_max, 10*w0, ymax))),
						"Lorentz":					(lambda *x: lineshape("Lorentz", 0, *x[:-1])+x[-1], (x0, y0, w0, ymean),      ((xmin, amp_min, 0, ymin),    (xmax, amp_max, 10*w0, ymax))),
						"Voigt":					(lambda *x: lineshape("Voigt", 0, *x[:-1])+x[-1],   (x0, y0, w0, w0, ymean),  ((xmin, amp_min, 0, 0, ymin), (xmax, amp_max, 5*w0, 5*w0, ymax))),
						"Gauss 1st Derivative":		(lambda *x: lineshape("Gauss", 1, *x[:-1])+x[-1],   (x0, y0, w0, ymean),      ((xmin, amp_min, 0, ymin),    (xmax, amp_max, 10*w0, ymax))),
						"Lorentz 1st Derivative":	(lambda *x: lineshape("Lorentz", 1, *x[:-1])+x[-1], (x0, y0, w0, ymean),      ((xmin, amp_min, 0, ymin),    (xmax, amp_max, 10*w0, ymax))),
						"Voigt 1st Derivative":		(lambda *x: lineshape("Voigt", 1, *x[:-1])+x[-1],   (x0, y0, w0, w0, ymean),  ((xmin, amp_min, 0, 0, ymin), (xmax, amp_max, 5*w0, 5*w0, ymax))),
						"Gauss 2nd Derivative":		(lambda *x: lineshape("Gauss", 2, *x[:-1])+x[-1],   (x0, y0, w0, ymean),      ((xmin, amp_min, 0, ymin),    (xmax, amp_max, 10*w0, ymax))),
						"Lorentz 2nd Derivative":	(lambda *x: lineshape("Lorentz", 2, *x[:-1])+x[-1], (x0, y0, w0, ymean),      ((xmin, amp_min, 0, ymin),    (xmax, amp_max, 10*w0, ymax))),
						"Voigt 2nd Derivative":		(lambda *x: lineshape("Voigt", 2, *x[:-1])+x[-1],   (x0, y0, w0, w0, ymean),  ((xmin, amp_min, 0, 0, ymin), (xmax, amp_max, 5*w0, 5*w0, ymax))),
					}.get(fitfunction)
				else:
					function, p0, bounds = {
						"Gauss":					(lambda *x: lineshape("Gauss", 0, *x),   (x0, y0, w0),      ((xmin, amp_min, 0),    (xmax, amp_max, 10*w0))),
						"Lorentz":					(lambda *x: lineshape("Lorentz", 0, *x), (x0, y0, w0),      ((xmin, amp_min, 0),    (xmax, amp_max, 10*w0))),
						"Voigt":					(lambda *x: lineshape("Voigt", 0, *x),   (x0, y0, w0, w0),  ((xmin, amp_min, 0, 0), (xmax, amp_max, 5*w0, 5*w0))),
						"Gauss 1st Derivative":		(lambda *x: lineshape("Gauss", 1, *x),   (x0, y0, w0),      ((xmin, amp_min, 0),    (xmax, amp_max, 10*w0))),
						"Lorentz 1st Derivative":	(lambda *x: lineshape("Lorentz", 1, *x), (x0, y0, w0),      ((xmin, amp_min, 0),    (xmax, amp_max, 10*w0))),
						"Voigt 1st Derivative":		(lambda *x: lineshape("Voigt", 1, *x),   (x0, y0, w0, w0),  ((xmin, amp_min, 0, 0), (xmax, amp_max, 5*w0, 5*w0))),
						"Gauss 2nd Derivative":		(lambda *x: lineshape("Gauss", 2, *x),   (x0, y0, w0),      ((xmin, amp_min, 0),    (xmax, amp_max, 10*w0))),
						"Lorentz 2nd Derivative":	(lambda *x: lineshape("Lorentz", 2, *x), (x0, y0, w0),      ((xmin, amp_min, 0),    (xmax, amp_max, 10*w0))),
						"Voigt 2nd Derivative":		(lambda *x: lineshape("Voigt", 2, *x),   (x0, y0, w0, w0),  ((xmin, amp_min, 0, 0), (xmax, amp_max, 5*w0, 5*w0))),
					}.get(fitfunction)


				try:
					popt, pcov = optimize.curve_fit(function, exp_xs, exp_ys, p0=p0, bounds=bounds)
				except Exception as E:
					popt, pcov = optimize.curve_fit(function, exp_xs, exp_ys, p0=p0, bounds=bounds)
				perr = np.sqrt(np.diag(pcov))
				fit_ys = function(fit_xs, *popt)

				xmiddle = popt[0]
				xuncert = perr[0]

			if index:
				ax = self.axs["ax"][index]
				self.fitcurve = ax.plot(fit_xs, fit_ys, color=main.config["color_fit"], alpha=0.7, linewidth=1)[0]
				self.fitline = ax.axvline(x=xmiddle, color=main.config["color_fit"], ls="--", alpha=1, linewidth=1)

		except Exception as E:
			self.fitcurve = None
			self.fitline = None
			main.notification(f"<span style='color:#eda711;'>WARNING</span>: The fitting failed with the following error message : {str(E)}")
			raise

		return(xmiddle, xuncert)

	def assign(self, index=None, init_values={}):
		lin_dict = {
			"qnu1":			pyckett.SENTINEL,
			"qnu2":			pyckett.SENTINEL,
			"qnu3":			pyckett.SENTINEL,
			"qnu4":			pyckett.SENTINEL,
			"qnu5":			pyckett.SENTINEL,
			"qnu6":			pyckett.SENTINEL,
			"qnl1":			pyckett.SENTINEL,
			"qnl2":			pyckett.SENTINEL,
			"qnl3":			pyckett.SENTINEL,
			"qnl4":			pyckett.SENTINEL,
			"qnl5":			pyckett.SENTINEL,
			"qnl6":			pyckett.SENTINEL,
			"x":			pyckett.SENTINEL,
			"error":		0,
			"weight":		1,
			"comment":		main.config["fit_comment"],
			"filename":		"__lin__",
			"xpre":			None,
		}

		lin_dict.update(init_values)

		error_value = main.config["fit_error"]

		if error_value > 0:
			lin_dict["error"] = error_value
		elif error_value >= -1:
			lin_dict["error"] = abs(lin_dict["x"] - lin_dict["xpre"])
		elif error_value >= -2:
			resp, rc = QInputDialog.getText(self, 'Set error', 'Error:')
			if rc:
				try:
					lin_dict["error"] = float(resp)
				except ValueError:
					main.notification("<span style='color:#eda711;'>WARNING</span>: Did not understand the given value for the error. The line was not assigned.")
					return
			else:
				return
		elif error_value >= -3:
			# This is the case, where the error is taken from the fitfunction (provided via init_values)
			pass

		if SeriesFitWindow in main.open_windows and main.open_windows[SeriesFitWindow].isVisible() and main.config["seriesfitwindow_greedy"]:
			lin_dict = {key: value for key, value in lin_dict.items() if key in main.ser_df.columns}
			lin_df = main.ser_df
			assigned_marker = False
		else:
			lin_dict = {key: value for key, value in lin_dict.items() if key in main.lin_df.columns}
			lin_df = main.new_df
			assigned_marker = True

		lin_df.reset_index(drop=True, inplace=True)
		new_ind = len(lin_df.index)
		lin_df.loc[new_ind] = lin_dict
		main.signalclass.assignment.emit()
		self.lastassignment = None

		if index and assigned_marker:
			ax = self.axs["ax"][index]

			annotation = self.axs["annotation"][index]
			reference = self.get_series_reference(index[1])
			if annotation and reference["method"] == "Transition":
				annotation.set_color(main.config["color_lin"])

			lin_plot = self.axs["lin_plot"][index]
			if lin_plot:
				offsets = lin_plot.get_offsets()
				offsets = np.concatenate([offsets, np.array((lin_dict["x"], 0), ndmin=2)])
				lin_plot.set_offsets(offsets)

				colors = lin_plot.get_facecolor()
				color = matplotlib.colors.to_rgba(main.config["color_lin"])
				colors = np.concatenate([colors, np.array((color, ), ndmin=2)])
				lin_plot.set_color(colors)

			with locks["axs"]:
				main.signalclass.drawplot.emit()

	@synchronized_d(locks["lin_df"])
	@synchronized_d(locks["ser_df"])
	def manual_assign(self):
		if self.lastassignment is None:
			main.notification("No saved assignment.")
			return
		index, dict_ = self.lastassignment
		
		if self.check_blends(index, dict_):
			pass
		else:
			self.assign(index, dict_)
		self.lastassignment = None


	def get_current_plot(self):
		lcp = self.lcp
		shape = self.axs["ax"].shape
		if self.is_valid_plot(lcp):
			return(lcp)
		else:
			return((0, 0))

	def is_valid_plot(self, i):
		shape = self.axs["ax"].shape
		if 0 <= i[0] < shape[0] and 0 <= i[1] < shape[1]:
			return(True)
		else:
			return(False)
			
	def create_plots(self):
		thread = threading.Thread(target=self.create_plots_core)
		with locks["currThread"]:
			thread.start()
			self.create_plots_id = thread.ident
		return(thread)

	@working_d
	@synchronized_d(locks["axs"])
	def create_plots_core(self):
		with locks["currThread"]:
			ownid = threading.current_thread().ident

		try:
			rows = max(main.config["plot_rows"], 1)
			cols = max(main.config["plot_cols"], 1)

			# The following performed better than self.fig.clf()
			for ax in self.fig.get_axes():
				self.fig.delaxes(ax)

			earlyreturn(ownid, self.create_plots_id)

			tmp = self.fig.subplots(rows, cols, gridspec_kw=main.config["plot_matplotlibkwargs"], squeeze=False)

			earlyreturn(ownid, self.create_plots_id)

			self.axs = {
				"ax": np.empty(tmp.shape, dtype=object),
				"exp_plot": np.empty(tmp.shape, dtype=object),
				"cat_plot": np.empty(tmp.shape, dtype=object),
				"lin_plot": np.empty(tmp.shape, dtype=object),
				"span": np.empty(tmp.shape, dtype=object),
				"annotation": np.empty(tmp.shape, dtype=object),
				"offset": np.zeros(tmp.shape),
				"width": np.full(tmp.shape, main.config["plot_width"], dtype=np.float64),
			}

			earlyreturn(ownid, self.create_plots_id)

			for i, row in enumerate(tmp):
				for j, ax in enumerate(row):
					exp_coll = matplotlib.collections.LineCollection(np.zeros(shape=(0,2,2)), colors=main.config["color_exp"])
					cat_coll = matplotlib.collections.LineCollection(np.zeros(shape=(0,2,2)), colors=main.config["color_cat"])
					span = matplotlib.widgets.SpanSelector(ax, lambda xmin, xmax, i=i, j=j:self.on_range(xmin, xmax, (i, j)), 'horizontal', useblit=True)

					self.axs["ax"][i, j] = ax
					self.axs["exp_plot"][i, j] = exp_coll
					self.axs["cat_plot"][i, j] = cat_coll
					self.axs["lin_plot"][i, j] = ax.scatter([], [], color=main.config["color_lin"], marker="*", zorder=100)
					self.axs["span"][i, j] = span

					ax.add_collection(exp_coll)
					ax.add_collection(cat_coll)

					ax.yaxis.set_visible(False)
					ax.xaxis.set_visible(False)

			earlyreturn(ownid, self.create_plots_id)

			for bottomax in self.axs["ax"][-1, :]:
				bottomax.xaxis.set_visible(True)

			earlyreturn(ownid, self.create_plots_id)

			main.signalclass.createdplots.emit()
			self.set_data()
		except CustomError as E:
			pass

	def manual_draw(self):
		self.set_data().join()
		with locks["axs"]:
			main.signalclass.drawplot.emit()

	def set_data(self):
		thread = threading.Thread(target=self.set_data_core)
		with locks["currThread"]:
			thread.start()
			self.set_data_id = thread.ident
		return(thread)

	@working_d
	@synchronized_d(locks["axs"])
	def set_data_core(self):
		with locks["currThread"]:
			ownid = threading.current_thread().ident

		try:
			if not main.config["flag_automatic_draw"]:
				return

			earlyreturn(ownid, self.set_data_id)

			# set x-ranges
			cat_df = main.get_visible_data("cat", scale=False)
			xpos, qns = self.get_positions(return_qns=True, cat_df=cat_df)
			widths = self.get_widths()
			widths += (widths == 0)*10

			offsets = self.get_offset(all=True)
			xmax = xpos + offsets + widths/2
			xmin = xpos + offsets - widths/2

			earlyreturn(ownid, self.set_data_id)

			# set ticks for bottom row
			i = -1
			for j in range(self.axs["ax"].shape[1]):
				ax = self.axs["ax"][i, j]
				ticks = np.linspace(xmin[i, j], xmax[i, j], main.config["plot_ticks"])

				reference = self.get_series_reference(j)
				if reference["method"] == "Fortrat":
					ticklabels = symmetric_ticklabels(ticks/2/reference["fortrat"]["jstart"])
				elif main.config["plot_offsetticks"] == 1:
					ticklabels = symmetric_ticklabels(ticks-xpos[i, j])
				elif main.config["plot_offsetticks"] == 0:
					ticklabels = symmetric_ticklabels(ticks)
				else:
					ticklabels = [f"{x:.2e}".replace("e+00", "").rstrip("0").rstrip(".") for x in ticks]

				if j!=0 and len(ticks)>1:
					ticks = ticks[1:]
					ticklabels = ticklabels[1:]
				ax.set_xticks(ticks)
				ax.set_xticklabels(ticklabels)

			earlyreturn(ownid, self.set_data_id)

			# set data and set y-ranges of data
			bins = main.config["plot_bins"]
			nobinning = main.config["plot_skipbinning"]
			scaling = main.config["plot_yscale"]

			datatypes = ("exp", "cat", "lin")
			dataframes = [main.get_visible_data("exp", scale=False), cat_df, main.get_visible_data("lin", scale=False)]
			files_dicts = [main.config[f"files_{type}"] for type in datatypes]

			minindices = {type: dataframe["x"].searchsorted(xmin, side="left") for type, dataframe in zip(datatypes, dataframes)}
			maxindices = {type: dataframe["x"].searchsorted(xmax, side="right") for type, dataframe in zip(datatypes, dataframes)}

			scalingfactordicts = {type: {file: main.config[f"files_{type}"][file].get("scale", 1) for file in main.config[f"files_{type}"].keys()} for type in ("exp", "cat")}

			for i in range(self.axs["ax"].shape[0]):
				for j in range(self.axs["ax"].shape[1]):
					earlyreturn(ownid, self.set_data_id)
					ax = self.axs["ax"][i, j]
					ax.set_xlim(xmin[i,j], xmax[i,j])

					for datatype, dataframe, files in zip(datatypes, dataframes, files_dicts):
						minindex, maxindex = minindices[datatype][i][j], maxindices[datatype][i][j]
						dataframe = dataframe.iloc[minindex:maxindex].copy()

						# Tested to thread this part as the bin_data function uses pandas functions that release the Global Interpreter Lock (GIL)
						# However, threading would only help in engineered cases where the data in each plot is huge (meaning high binning parameter but
						# even higher number of points)
						# It was decided, that the added complexity was not worth the benefit for these extreme cases. Additionally, the overhead of threads
						# would reduce performance for the majority of use cases
						binwidth = (xmax[i, j]-xmin[i, j]) / bins
						if len(dataframe) > max(bins, nobinning)  and binwidth != 0:
							dataframe = bin_data(dataframe, binwidth, (xmin[i, j], xmax[i, j]))

						xs = dataframe["x"].to_numpy()
						if datatype == "lin":
							ys = 0*xs
						else:
							ys = dataframe["y"]*dataframe["filename"].replace(scalingfactordicts[datatype]).to_numpy()

						if datatype == "exp":
							if scaling == "Per Plot":
								if len(dataframe):
									yrange_exp = [ys.min(), ys.max()]
								else:
									yrange_exp = [-1, 1]
							if main.config["plot_expasstickspectrum"]:
								segs = np.array(((xs, xs), (ys*0, ys))).T
								colors = create_colors(dataframe, files)
							else:
								# @Luis: optimize this
								filenames = dataframe["filename"].to_numpy()
								unique_filenames = np.unique(filenames)

								segs = []
								colors = []
								for unique_filename in unique_filenames:
									mask = (filenames == unique_filename)
									tmp_xs, tmp_ys = xs[mask], ys[mask]

									segs.append(np.array(((tmp_xs[:-1], tmp_xs[1:]), (tmp_ys[:-1], tmp_ys[1:]))).T)
									colors.extend([files[unique_filename]["color"]]*sum(mask))

								if segs:
									segs = np.concatenate(segs)

							self.axs["exp_plot"][i, j].set(segments=segs, colors=colors)
						elif datatype == "cat":
							if scaling == "Per Plot":
								if len(dataframe):
									yrange_cat = [ys.min(), ys.max()]
								else:
									yrange_cat = [-1, 1]
								ys = ys*yrange_exp[1]/yrange_cat[1]
							elif scaling in ["Global", "Custom"]:
								ys = ys*main.config["plot_expcat_factor"]*10**main.config["plot_expcat_exponent"]
							segs = np.array(((xs, xs), (ys*0, ys))).T
							# segs = (((xs[i], 0),(xs[i], ys[i])) for i in range(len(xs)))
							colors = create_colors(dataframe, files, xpos[i, j])

							self.axs["cat_plot"][i, j].set(segments=segs, colors=colors)
						elif datatype == "lin":
							tuples = list(zip(xs,ys))
							tuples = tuples if len(tuples)!=0 else [[None,None]]
							colors = create_colors(dataframe, files, lin=True)

							self.axs["lin_plot"][i, j].set_offsets(tuples)
							self.axs["lin_plot"][i, j].set_color(colors)

					if scaling == "Per Plot":
						yrange = yrange_exp
					elif scaling == "Global":
						yrange = main.yrange_exp
					else:
						yrange = (main.config["plot_yscale_min"], main.config["plot_yscale_max"])

					margin = main.config["plot_ymargin"]

					yrange = [yrange[0]-margin*(yrange[1]-yrange[0]), yrange[1]+margin*(yrange[1]-yrange[0])]
					if np.isnan(yrange[0]) or np.isnan(yrange[1]) or yrange[0] == yrange[1]:
						yrange = [-1,+1]
					ax.set_ylim(yrange)


			earlyreturn(ownid, self.set_data_id)
			self.plot_annotations(xpos, qns)

			earlyreturn(ownid, self.set_data_id)
			main.signalclass.drawplot.emit()

		except CustomError as E:
			pass

	def wheelEvent(self, event):
		steps = event.angleDelta().y() // 120
		self.set_width(2 ** -steps, absolute=False)

	def change_fitcolor(self):
		color = QColorDialog.getColor(initial=QColor(rgbt_to_trgb(main.config["color_fit"])), options=QColorDialog.ColorDialogOption.ShowAlphaChannel)
		color = Color(trgb_to_rgbt(color.name(QColor.NameFormat.HexArgb)))
		main.config["color_fit"] = color

	def plot_number(self):
		resp, rc = QInputDialog.getText(self, 'How many plots do you want: ', 'Number:')
		if not rc:
			return
		resp = resp.split("x")
		try:
			if len(resp) == 1:
				main.config["plot_rows"] = int(resp[0])
			elif len(resp) == 2:
				main.config["plot_rows"] = int(resp[0])
				main.config["plot_cols"] = int(resp[1])
		except ValueError:
			main.notification("<span style='color:#eda711;'>WARNING</span>: The entered value was not understood.")
			return

	def offset_dialog(self):
		resp, rc = QInputDialog.getText(self, 'Go to frequency', 'Frequency:')
		if not rc:
			return

		try:
			xnew = float(resp)
		except ValueError:
			main.notification("<span style='color:#eda711;'>WARNING</span>: The entered value could not be interpreted as a number.")
			return

		absolute = not main.config["plot_relativegoto"]
		self.set_offset(xnew, absolute=absolute)

	def width_dialog(self):
		resp, rc = QInputDialog.getText(self, 'Set view width', 'Width:')
		if not rc:
			return
		try:
			self.set_width(float(resp))
		except ValueError:
			main.notification("<span style='color:#eda711;'>WARNING</span>: The entered value could not be interpreted as a number.")
			return




##
## Enhanced DockWidget Class
##
class EQDockWidget(QDockWidget):
	def __init__(self, parent=None):
		super().__init__(main.mainwindow)
		self.setObjectName(self.__class__.__name__)

		main.mainwindow.addDockWidget(Qt.DockWidgetArea(2), self)
		QShortcut("Esc", self).activated.connect(self.close)

class ConfigurePlotsWindow(EQDockWidget):
	def __init__(self, parent=None):
		super().__init__(parent)
		self.setWindowTitle("Configure Plots")

		mainwidget = QGroupBox()
		layout = QVBoxLayout()
		self.setWidget(mainwidget)
		mainwidget.setLayout(layout)

		hbox1 = QHBoxLayout()
		binninglabel = QQ(QLabel, text="Binning: ")
		binningfield = QQ(QSpinBox, "plot_bins", range=(1, None), minWidth=120)
		[hbox1.addWidget(binninglabel), hbox1.addWidget(binningfield), hbox1.addStretch(1)]

		checkbox_coupled  = QQ(QCheckBox, "plot_coupled", text="Plots are coupled")
		checkbox_annotate = QQ(QCheckBox, "plot_annotate", text="Annotate plots")


		[layout.addLayout(hbox1), layout.addItem(QSpacerItem(5, 5, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding))]
		[layout.addWidget(widget) for widget in (checkbox_coupled, checkbox_annotate, )]

		checkbox_scale = QQ(QComboBox, "plot_yscale", options=("Per Plot", "Global", "Custom"), change=lambda x: main.signalclass.updateplot.emit())
		spinbox_scalemin = QQ(QDoubleSpinBox, "plot_yscale_min", range=(None, None), minWidth=120, change=lambda x: main.signalclass.updateplot.emit())
		spinbox_scalemax = QQ(QDoubleSpinBox, "plot_yscale_max", range=(None, None), minWidth=120, change=lambda x: main.signalclass.updateplot.emit())
		spinbox_scalecatfac = QQ(QDoubleSpinBox, "plot_expcat_factor", range=(0, None), minWidth=120, change=lambda x: main.signalclass.updateplot.emit())
		spinbox_scalecatexp = QQ(QSpinBox, "plot_expcat_exponent", range=(None, None), prefix="*10^", minWidth=120, change=lambda x: main.signalclass.updateplot.emit())

		scaleMinLabel = QLabel("y-min:  ")
		scaleMaxLabel = QLabel("y-max:  ")
		scaleCatLabel = QLabel("y-Exp/y-Cat:  ")

		checkbox_scale.currentTextChanged.connect(lambda x: [
			spinbox_scalecatfac.setVisible(x!="Per Plot"),
			spinbox_scalecatexp.setVisible(x!="Per Plot"),
			spinbox_scalemin.setVisible(x=="Custom"),
			spinbox_scalemax.setVisible(x=="Custom"),
			scaleMaxLabel.setVisible(x=="Custom"),
			scaleMinLabel.setVisible(x=="Custom"),
			scaleCatLabel.setVisible(x!="Per Plot"),
		])
		checkbox_scale.currentTextChanged.emit(main.config["plot_yscale"])

		grid = QGridLayout()
		grid.addWidget(QLabel("y-scale:  "), 0, 0)
		grid.addWidget(checkbox_scale, 0, 1)
		grid.addWidget(scaleCatLabel, 1, 0)
		grid.addWidget(spinbox_scalecatfac, 1, 1)
		grid.addWidget(spinbox_scalecatexp, 1, 2, 1, 2)
		grid.addWidget(scaleMinLabel, 2, 0)
		grid.addWidget(spinbox_scalemin, 2, 1)
		grid.addWidget(scaleMaxLabel, 2, 2)
		grid.addWidget(spinbox_scalemax, 2, 3)
		grid.setColumnStretch(7, 10)
		grid.setRowStretch(2, 10)
		layout.addItem(QSpacerItem(5, 5, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding))
		layout.addLayout(grid)

		layout.addStretch(1)

class ReferenceSeriesWindow(EQDockWidget):
	def __init__(self, parent=None):
		super().__init__(parent)
		self.setWindowTitle("Reference Series")

		mainwidget = QGroupBox()
		layout = QVBoxLayout()
		self.setWidget(mainwidget)
		mainwidget.setLayout(layout)

		self.tabs = QTabWidget()
		layout.addWidget(self.tabs)

		self.tabs.setTabsClosable(True)
		self.tabs.setMovable(True)
		self.tabs.setDocumentMode(True)

		self.tabs.setCornerWidget(QQ(QToolButton, text="Dupl.", tooltip="Duplicate current tab", change=self.duplicate_tab), Qt.Corner.TopRightCorner)
		self.tabs.tabCloseRequested.connect(self.close_tab)
		self.tabs.tabBarDoubleClicked.connect(self.renameoradd_tab)
		self.tabs.setCurrentIndex(main.config["series_currenttab"])
		self.tabs.currentChanged.connect(lambda x: main.config.__setitem__("series_currenttab", x))
		main.signalclass.createdplots.connect(self.min_series)

		self.set_state(main.config["series_references"])
		if not self.tabs.count():
			self.add_tab()

		self.tabs.currentChanged.connect(lambda x: self.changed())

	def close_tab(self, index):
		if self.tabs.count() == 1:
			return
		tab = self.tabs.widget(index)
		tab.deleteLater()
		self.tabs.removeTab(index)
		self.get_state()

	def renameoradd_tab(self, index):
		if index == -1:
			self.add_tab()
		elif self.tabs.widget(index) != 0:
			text, ok = QInputDialog().getText(self, "Tab Name","Enter the Tabs Name:")
			if ok and text:
				self.tabs.setTabText(index, text)
				self.tabs.widget(index).state["title"] = text

	def add_tab(self, init_values={}):
		title = init_values.get("title", "Series")
		tmp = ReferenceSelector(init_values, self)
		self.tabs.addTab(tmp, title)
		self.changed()

	def duplicate_tab(self):
		values = self.tabs.widget(self.tabs.currentIndex()).get_values()
		self.add_tab(values)

	def min_series(self):
		numberoftabs = self.tabs.count()
		locked = main.config["flag_referencenumberlocked"]
		diff = main.config["plot_cols"] - numberoftabs

		if diff > 0:
			for i in range(diff):
				self.duplicate_tab()
		elif diff < 0 and locked:
			for i in range(-diff):
				self.close_tab(numberoftabs - i -1)

	def set_state(self, values):
		for i in range(self.tabs.count()):
			tab = self.tabs.widget(i)
			tab.deleteLater()
		self.tabs.clear()
		
		for tabdata in values:
			self.add_tab(tabdata)

	def get_state(self):
		values = []
		for i in range(self.tabs.count()):
			tmp = self.tabs.widget(i).get_values()
			values.append(tmp)
		return(values)

	def changed(self):
		tmp = self.get_state()
		main.config["series_references"] = self.get_state()

class CatalogWindow(EQDockWidget):
	def __init__(self, parent=None):
		super().__init__(parent)
		self.setWindowTitle("Catalog of newly assigned Lines")

		mainwidget = QGroupBox()
		layout = QVBoxLayout()
		self.setWidget(mainwidget)
		mainwidget.setLayout(layout)

		button_savecatalog = QQ(QToolButton, text="Save", change=lambda: main.save_lines_lin())
		button_delete = QQ(QToolButton, text="X", change=main.catalog_table_delete)
		button_deleteall = QQ(QToolButton, text="Del All", change=lambda: main.catalog_table_delete(all=True))
		button_addrow = QQ(QToolButton, text="+", change=lambda x: addemptyrow_inplace(main.new_df, self.catalogModel))
		button_tight = QQ(QToolButton, text="Tight", change=lambda x: self.catalogTable.resizeColumnsToContents())
		checkbox_appendonsave = QQ(QCheckBox, "flag_appendonsave", text="Append", tooltip="Append to file if checked or overwrite content if unchecked")
		errorlabel = QQ(QLabel, text="Default Uncertainty: ", tooltip="Positive Values are absolute Values, -1 for obs-calc, -2 for dialog, -3 for StdDev from Fit")
		errorfield = QQ(QDoubleSpinBox, "fit_error", range=(-3, None), minWidth=120, singlestep=main.config["fit_uncertaintystep"])

		headers = ["U1", "U2", "U3", "U4", "U5", "U6", "L1", "L2", "L3", "L4", "L5", "L6", "Freq", "Unc.", "Weight", "Comment"]
		str_columns = [headers.index("Comment")]

		self.catalogTable = QTableView()
		self.catalogModel = CustomTableModel(main.new_df, headers, ["filename"])
		self.catalogTable.setModel(self.catalogModel)
		self.catalogTable.resizeColumnsToContents()

		main.signalclass.assignment.connect(self.scroll_bottom)

		main.config.register("series_qns", self.update_columns_visibility)
		self.update_columns_visibility()

		buttonsBox = QHBoxLayout()
		[buttonsBox.addWidget(button_delete), buttonsBox.addWidget(button_deleteall), buttonsBox.addWidget(button_addrow),
		 buttonsBox.addStretch(2), buttonsBox.addWidget(checkbox_appendonsave), buttonsBox.addWidget(button_savecatalog),
		 buttonsBox.addWidget(button_tight)]
		layout.addLayout(buttonsBox)
		layout.addWidget(self.catalogTable)
		buttonsBox = QHBoxLayout()
		[buttonsBox.addWidget(errorlabel), buttonsBox.addStretch(1), buttonsBox.addWidget(errorfield)]
		layout.addLayout(buttonsBox)

	def scroll_bottom(self):
		self.catalogTable.selectRow(len(self.catalogModel.data)-1)
		self.catalogModel.update()
		self.catalogTable.scrollToBottom()

	def update_columns_visibility(self):
		qns = main.config["series_qns"]
		for i in range(6):
			self.catalogTable.setColumnHidden(i,   i>=qns)
			self.catalogTable.setColumnHidden(i+6, i>=qns)

class LogWindow(EQDockWidget):
	def __init__(self, parent=None):
		super().__init__(parent)
		self.setWindowTitle("Log")

		mainwidget = QGroupBox()
		layout = QVBoxLayout()
		self.setWidget(mainwidget)
		mainwidget.setLayout(layout)

		self.log_area = QTextEdit()
		self.log_area.setReadOnly(True)
		self.log_area.setMinimumHeight(50)

		main.signalclass.writelog.connect(lambda text: self.writelog(text))
		layout.addWidget(self.log_area)

	def writelog(self, text):
		if main.config["flag_alwaysshowlog"]:
			self.setVisible(True)
			self.raise_()

		tmp = self.log_area.toPlainText()
		tmp = tmp.split("\n")
		if len(tmp)-1 > main.config["flag_logmaxrows"]:
			self.log_area.setText("\n".join(tmp[-main.config["flag_logmaxrows"]:]))

		self.log_area.append(text)
		sb = self.log_area.verticalScrollBar()
		sb.setValue(sb.maximum())

class HoverWindow(EQDockWidget):
	def __init__(self, parent=None):
		super().__init__(parent)
		self.setWindowTitle("Hover")

		mainwidget = QGroupBox()
		layout = QVBoxLayout()
		self.setWidget(mainwidget)
		mainwidget.setLayout(layout)

		self.log_area = QTextEdit()
		self.log_area.setReadOnly(True)
		self.log_area.setMinimumHeight(50)

		main.signalclass.writehover.connect(lambda text: self.log_area.setText(text))
		layout.addWidget(self.log_area)

class QuoteWindow(EQDockWidget):
	def __init__(self, parent=None):
		super().__init__(parent)
		self.setWindowTitle("Quote")

		mainwidget = QGroupBox()
		layout = QVBoxLayout()
		self.setWidget(mainwidget)
		mainwidget.setLayout(layout)

		self.quote = QQ(QLabel, wordwrap=True, align=Qt.AlignmentFlag.AlignCenter)
		self.new_quote()
		layout.addWidget(self.quote)
		self.setVisible(False)

	def toggle(self):
		self.new_quote()
		self.setVisible(not self.isVisible())

	def new_quote(self):
		quotes = json.loads(quotes_str)
		quote = quotes[random.randint(0,len(quotes)-1)]
		self.quote.setText(quote)
		self.quote.setSizePolicy(QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Minimum)




##
## Enhanced Window Class
##
class EQWidget(QWidget):
	def __init__(self, id, parent=None):
		self.id = id
		geometry = main.config.get(f"windowgeometry_{self.id}")
		super().__init__(parent)

		if geometry:
			if isinstance(geometry, str):
				geometry = json.loads(geometry)
			self.setGeometry(*geometry)

		QShortcut("Esc", self).activated.connect(self.close)

	@synchronized_d(locks["windows"])
	def closeEvent(self, *args, **kwargs):
		main.config[f"windowgeometry_{self.id}"] = self.geometry().getRect()
		return super().closeEvent(*args, **kwargs)

	def resizeEvent(self, event):
		main.config[f"windowgeometry_{self.id}"] = self.geometry().getRect()

	def moveEvent(self, event):
		main.config[f"windowgeometry_{self.id}"] = self.geometry().getRect()

	def show(self, *args, **kwargs):
		screen_box = self.screen().geometry()
		widget_top_left = self.geometry().topLeft()
		widget_bottom_right = self.geometry().bottomRight()
		
		if not (screen_box.contains(widget_top_left) and screen_box.contains(widget_bottom_right)):
			primary_screen = QApplication.instance().primaryScreen()
			self.move(primary_screen.geometry().center()- self.rect().center())
		
		return(super().show(*args, **kwargs))

class CreditsWindow(EQWidget):
	def __init__(self, id, parent=None):
		super().__init__(id, parent)
		self.setWindowTitle("Credits")

		global CREDITSSTRING
		layout = QVBoxLayout()
		layout.addWidget(QQ(QLabel, text=CREDITSSTRING, align=Qt.AlignmentFlag.AlignCenter, wordwrap=True, minHeight=300, minWidth=500))
		self.setLayout(layout)

class FileWindow(EQWidget):
	def __init__(self, id, parent=None):
		super().__init__(id, parent)
		self.setWindowTitle("Files Window")

		self.setAcceptDrops(True)
		self.tabs = QTabWidget()
		tmplayout = QVBoxLayout()
		tmplayout.addWidget(self.tabs)
		self.setLayout(tmplayout)

		keys = ("exp", "cat", "lin")
		self.widgets = {key: {} for key in keys}
		self.layouts = {key: self.create_layout(key, initial=True) for key in keys}
		for label, layout in self.layouts.items():
			tmpwidget = QWidget()
			tmpwidget.setLayout(layout)
			self.tabs.addTab(tmpwidget, label.capitalize())

		main.signalclass.fileschanged.connect(self.update)

	def dragEnterEvent(self, event):
		if event.mimeData().hasUrls():
			event.accept()
		else:
			event.ignore()

	def dropEvent(self, event):
		main.mainwindow.dropEvent(event)

	def update(self, type=None):
		if type is None:
			types = ("exp", "cat", "lin")
		else:
			types = (type, )


		for type in types:
			filesgrid = self.widgets[f"{type}_filesgrid"]

			scrollarea = self.widgets.get(f"{type}_scrollarea")
			if scrollarea:
				tmp = (scrollarea.verticalScrollBar().value(), scrollarea.horizontalScrollBar().value())
			else:
				tmp = (0, 0)

			# Delete existing widgets
			for key, value in self.widgets[type].items():
				if not (key.startswith("__") and key.endswith("__")):
					for widget in value.values():
						widget.deleteLater()

			self.widgets[type] = {key: value for key, value in self.widgets[type].items() if (key.startswith("__") and key.endswith("__"))}

			if type == "exp":
				actions = ("label", "colorinput", "colorpicker", "scale", "hide", "delete", "reread")
			elif type == "cat":
				actions = ("label", "colorinput", "colorpicker", "scale", "hide", "delete", "reread")
			elif type == "lin":
				actions = ("label", "colorinput", "colorpicker", "hide", "delete", "reread")

			row_id = 0
			files = main.config[f"files_{type}"]

			for file in files:
				self.add_row(filesgrid, type, file, actions, row_id)
				row_id += 1

			filesgrid.setRowStretch(row_id, 1)

			scrollarea.verticalScrollBar().setValue(tmp[0])
			scrollarea.horizontalScrollBar().setValue(tmp[1])


	def create_layout(self, type, initial=False):
		if initial:
			layout = QVBoxLayout()

			buttonsbox = QHBoxLayout()
			scrollarea = QScrollArea()
			widget = QWidget()

			layout.addLayout(buttonsbox)

			buttonsbox.addWidget(QQ(QToolButton, text="Load", change=lambda x, type=type: main.load_file(type, add_files=True, keep_old=False)))
			buttonsbox.addWidget(QQ(QToolButton, text="Add", change=lambda x, type=type: main.load_file(type, add_files=True, keep_old=True)))
			buttonsbox.addWidget(QQ(QToolButton, text="Reread All", change=lambda x, type=type: main.load_file(type, reread=True, do_QNs=False)))
			buttonsbox.addWidget(QQ(QToolButton, text="Reset All", change=lambda x, type=type: self.reset_all(type)))
			buttonsbox.addWidget(QQ(QToolButton, text="Delete All", change=lambda x, type=type: self.delete_file(type)))
			buttonsbox.addStretch(1)

			filesgrid = QGridLayout()

			scrollarea.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
			scrollarea.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
			scrollarea.setWidgetResizable(True)
			scrollarea.setWidget(widget)
			widget.setLayout(filesgrid)

			self.widgets[f"{type}_filesgrid"] = filesgrid
			self.widgets[f"{type}_scrollarea"] = scrollarea


			if type == "exp":
				color = main.config.get('color_exp')
				file = "__exp__"

				topbox = QHBoxLayout()
				layout.addLayout(topbox)

				rowdict = {
					"label":			QQ(QLabel, text="Initial Color"),
					"colorinput":		QQ(QLineEdit, text=color, maxWidth=200, change=lambda x, file=file, type=type: self.change_color(type, file, inp=True)),
					"colorpicker":		QQ(QToolButton, text="CP", change=lambda x, file=file, type=type: self.change_color(type, file), stylesheet=f"background-color: {rgbt_to_trgb(color)}"),
				}
				for label, widget in rowdict.items():
					topbox.addWidget(widget)
				self.widgets[type][file] = rowdict

			elif type == "cat":
				labels = ("Default color", "Current transition")
				colors = (main.config.get('color_cat'), main.config.get('color_cur'))
				files = ("__cat__", "__cur__")

				topbox = QHBoxLayout()
				layout.addLayout(topbox)

				for label, color, file in zip(labels, colors, files):
					rowdict = {
						"label":			QQ(QLabel, text=label),
						"colorinput":		QQ(QLineEdit, text=color, maxWidth=200, change=lambda x, file=file, type=type: self.change_color(type, file, inp=True)),
						"colorpicker":		QQ(QToolButton, text="CP", change=lambda x, file=file, type=type: self.change_color(type, file), stylesheet=f"background-color: {rgbt_to_trgb(color)}"),
					}
					for label, widget in rowdict.items():
						topbox.addWidget(widget)
					self.widgets[type][file] = rowdict

			elif type == "lin":
				color = main.config.get('color_lin')
				file = "__lin__"
				hidden = main.config["flag_hidecatalog"]

				topbox = QHBoxLayout()
				layout.addLayout(topbox)

				rowdict = {
					"label":			QQ(QLabel, text="Assigned Markers"),
					"colorinput":		QQ(QLineEdit, text=color, maxWidth=200, change=lambda x, file=file, type=type: self.change_color(type, file, inp=True)),
					"colorpicker":		QQ(QToolButton, text="CP", change=lambda x, file=file, type=type: self.change_color(type, file), stylesheet=f"background-color: {rgbt_to_trgb(color)}"),
					"hide":				QQ(QToolButton, text="Show" if hidden else "Hide", change=lambda x, file=file, type=type: self.hide_file(type, file)),
				}
				for label, widget in rowdict.items():
					topbox.addWidget(widget)
				self.widgets[type][file] = rowdict

			layout.addWidget(scrollarea, 10)

		self.update(type)

		return(layout)

	def add_row(self, layout, type, file, actions, row_id):
		file_options = main.config[f"files_{type}"][file]
		color = file_options.get("color", "#ffffff")
		hidden = file_options.get("hidden")
		scale = file_options.get("scale", 1)

		rowdict = {
			"label":			QQ(QLabel, text=file, enabled=not hidden),
			"scale":			QQ(QDoubleSpinBox, value=scale, range=(None, None), change=lambda x, file=file, type=type: self.scale_file(type, file, x)),
			"colorinput":		QQ(QLineEdit, text=color, maxWidth=200, change=lambda x, file=file, type=type: self.change_color(type, file, inp=True)),
			"colorpicker":		QQ(QToolButton, text="CP", change=lambda x, file=file, type=type: self.change_color(type, file), stylesheet=f"background-color: {rgbt_to_trgb(color)}"),
			"hide":				QQ(QToolButton, text="Show" if hidden else "Hide", change=lambda x, type=type, file=file: self.hide_file(type, file)),
			"delete":			QQ(QToolButton, text="×", change=lambda x, type=type, file=file: self.delete_file(type, file), tooltip="Delete file"),
			"reread":			QQ(QToolButton, text="⟲", change=lambda x, type=type, file=file: main.load_file(type, add_files=[file], keep_old=True, do_QNs=False), tooltip="Reread File"),
		}

		for col_id, action in enumerate(actions):
			layout.addWidget(rowdict[action], row_id, col_id)
			layout.setRowStretch(row_id, 0)

		self.widgets[type][file] = rowdict

	def scale_file(self, type, file, scale):
		main.config[f"files_{type}"][file]["scale"] = scale
		main.plotwidget.set_data()

	def reset_all(self, type):
		files = main.config[f"files_{type}"]

		for file in files:
			if "scale" in files[file]:
				files[file]["scale"] = 1

			if "hidden" in files[file]:
				files[file]["hidden"] = False

			if "color" in files[file]:
				files[file]["color"] = main.config[f"color_{type}"]

		main.signalclass.fileschanged.emit()

	@working_d
	def delete_file(self, type, file=None):
		df = main.return_df(type)

		with locks[f"{type}_df"]:
			if file is None:
				main.config[f"files_{type}"].clear()
				df.drop(df.index, inplace=True)
			else:
				if file in main.config[f"files_{type}"]:
					del main.config[f"files_{type}"][file]
				df.drop(df[df["filename"]==file].index, inplace=True)

		main.load_file(type, keep_old=True, do_QNs=False)

	@synchronized_d(locks["axs"])
	def hide_file(self, type, file):
		if file == "__lin__":
			hidden = main.config["flag_hidecatalog"]
		else:
			hidden = main.config[f"files_{type}"][file].get("hidden", False)
		hidden = not hidden

		if file == "__lin__":
			main.config["flag_hidecatalog"] = hidden
		else:
			main.config[f"files_{type}"][file]["hidden"] = hidden

		if hidden:
			self.widgets[type][file]["label"].setEnabled(False)
			self.widgets[type][file]["hide"].setText("Show")
		else:
			self.widgets[type][file]["label"].setEnabled(True)
			self.widgets[type][file]["hide"].setText("Hide")

		main.signalclass.updateplot.emit()

	@synchronized_d(locks["axs"])
	def change_color(self, type, file, inp=False):
		color_input = self.widgets[type][file]["colorinput"].text()
		if inp:
			color = color_input
		else:
			color = QColorDialog.getColor(initial=QColor(rgbt_to_trgb(color_input)), options=QColorDialog.ColorDialogOption.ShowAlphaChannel)
			if color.isValid():
				color = trgb_to_rgbt(color.name(QColor.NameFormat.HexArgb))
			else:
				return

		try:
			color = Color(color)
		except CustomError:
			return

		self.widgets[type][file]["colorpicker"].setStyleSheet(f"background-color: {rgbt_to_trgb(color)}")
		if self.widgets[type][file]["colorinput"].text() != color:
			self.widgets[type][file]["colorinput"].setText(color)

		if file.startswith("__") and file.endswith("__"):
			tmp = file.strip("_")
			main.config[f"color_{tmp}"] = color
		else:
			main.config[f"files_{type}"][file]["color"] = color
		main.signalclass.updateplot.emit()

class LineshapeWindow(EQWidget):
	def __init__(self, id, parent=None):
		super().__init__(id, parent)
		self.setWindowTitle("Lineshape Window")

		class CustomPlotWidget(ProtPlot):
			def gui(self):
				super().gui()
				self.mod_line = self.ax.plot([], [], color = "blue")[0]

			def update_plot(self):
				cat_df = main.get_visible_data("cat", xrange=(self.center-self.width/2, self.center+self.width/2), binning=True)
				cat_xs = cat_df["x"].to_numpy()
				cat_ys = cat_df["y"].to_numpy()

				function = main.config["lineshapewindow_lineshape"]

				res_xs = np.linspace(self.center-self.width/2, self.center+self.width/2, main.config["lineshapewindow_xpoints"])
				res_ys = []

				args = [function, main.config["lineshapewindow_derivative"], res_xs, 0, 0]
				if function in ("Gauss", "Voigt"):
					args.append(main.config["lineshapewindow_gauss"])
				if function in ("Lorentz", "Voigt"):
					args.append(main.config["lineshapewindow_lorentz"])

				for x, y in zip(cat_xs, cat_ys):
					args[3] = x
					args[4] = y
					res_ys.append(lineshape(*args))

				if len(res_ys) > 0:
					res_ys = np.sum(np.array(res_ys), axis=0)
				else:
					res_ys = np.zeros(len(res_xs))
				res_ys = res_ys*main.config["lineshapewindow_scaling_factor"]*10**main.config["lineshapewindow_scaling_exponent"]

				self.mod_line.set_data(res_xs, res_ys)

				super().update_plot()

		layout = QVBoxLayout()
		self.setLayout(layout)
		self.plotWidget = CustomPlotWidget(parent=self)
		layout.addWidget(self.plotWidget, 1)

		tmp_layout = QGridLayout()
		layout.addLayout(tmp_layout)

		row_id = 0
		tmp_layout.addWidget(QQ(QLabel, text="Function: "), row_id, 0)
		tmp_layout.addWidget(QQ(QComboBox, "lineshapewindow_lineshape", minWidth=120, items=["Gauss", "Lorentz", "Voigt"]), row_id, 1)

		tmp_layout.addWidget(QQ(QLabel, text="Gauss: "), row_id, 2)
		tmp_layout.addWidget(QQ(QDoubleSpinBox, "lineshapewindow_gauss", minWidth=120), row_id, 3)

		row_id += 1
		tmp_layout.addWidget(QQ(QLabel, text="Derivative: "), row_id, 0)
		tmp_layout.addWidget(QQ(QSpinBox, "lineshapewindow_derivative", range=(0, 2), minWidth=120), row_id, 1)

		tmp_layout.addWidget(QQ(QLabel, text="Lorentz: "), row_id, 2)
		tmp_layout.addWidget(QQ(QDoubleSpinBox, "lineshapewindow_lorentz", minWidth=120), row_id, 3)

		row_id += 1
		tmp_layout.addWidget(QQ(QLabel, text="Scaling: "), row_id, 0)
		tmp2_layout = QHBoxLayout()
		tmp2_layout.addWidget(QQ(QDoubleSpinBox, "lineshapewindow_scaling_factor", range=(None, None)), 1)
		tmp2_layout.addWidget(QQ(QSpinBox, "lineshapewindow_scaling_exponent", prefix="*10^", range=(None, None)), 1)
		tmp_layout.addLayout(tmp2_layout, row_id, 1, 1, 3)

		tmp_layout.setColumnStretch(4, 1)

		main.config.register(("lineshapewindow_lineshape", "lineshapewindow_gauss", "lineshapewindow_lorentz", "lineshapewindow_derivative", "lineshapewindow_scaling_factor", "lineshapewindow_scaling_exponent"), self.plotWidget.update_plot)

	def activateWindow(self):
		self.plotWidget.from_current_plot()
		super().activateWindow()

class BlendedLinesWindow(EQWidget):
	def __init__(self, id, parent=None):
		super().__init__(id, parent)
		self.setWindowTitle("Blended Lines Window")

		self.peaks = []
		self.fit_values = None
		self.cid = None

		class CustomPlotWidget(ProtPlot):

			def gui(self):
				super().gui()
				self.fit_line = self.ax.plot([], [], color = main.config["blendedlineswindow_color_total"])[0]
				self.cat_line = matplotlib.collections.LineCollection(np.zeros(shape=(0,2,2)), color=main.config["color_cat"])
				self.ax.add_collection(self.cat_line)

				self.plot_parts = []

			@synchronized_d(locks["fitting"])
			def fit_peaks(self):
				with locks["currThread"]:
					ownid = threading.current_thread().ident

				try:
					main.signalclass.fitindicator.emit("<span style='color:#eda711;'>Working ...</span>")
					peaks = self.parent.peaks.copy()
					profile = main.config["blendedlineswindow_lineshape"]
					derivative = main.config["blendedlineswindow_derivative"]
					polynomrank = main.config["blendedlineswindow_polynom"]+1
					fixedwidth = main.config["blendedlineswindow_fixedwidth"]
					amplitude_direction = main.config["fit_amplitudedirection"]
					now = 2 if profile == "Voigt" else 1
					noa = 2 + now * (not fixedwidth)

					fitfunction = lambda x, *args, fixedwidth=fixedwidth: self.parent.fitfunction(x, profile, derivative, polynomrank, *args, fixedwidth=fixedwidth)
					fitfunction_withoutbaseline = lambda x, *args, fixedwidth=fixedwidth: self.parent.fitfunction(x, profile, derivative, 0, *args, fixedwidth=fixedwidth)

					for part in self.plot_parts:
						part.remove()
					self.plot_parts = []

					xrange = [self.center-self.width/2, self.center+self.width/2]

					earlyreturn(ownid, self.fit_thread_id)

					df_exp = main.get_visible_data("exp", xrange=xrange)
					exp_xs = df_exp["x"].to_numpy()
					exp_ys = df_exp["y"].to_numpy()

					earlyreturn(ownid, self.fit_thread_id)

					df_cat = main.get_visible_data("cat", xrange=xrange)
					if len(df_cat) != 0:
						cat_xs = df_cat["x"].to_numpy()
						cat_ys = df_cat["y"].to_numpy()
						if len(exp_ys):
							cat_ys = cat_ys*np.max(exp_ys)/np.max(cat_ys)

						segs = np.array(((cat_xs, cat_xs), (cat_ys*0, cat_ys))).T
						# segs = (((cat_xs[i],0),(cat_xs[i],cat_ys[i])) for i in range(len(cat_xs)))
						colors = create_colors(df_cat, main.config["files_cat"])
						self.cat_line.set(segments=segs, colors=colors)

					earlyreturn(ownid, self.fit_thread_id)

					xs = []
					ys = []
					ws = []
					exp_mean = 0
					if len(exp_ys) and (polynomrank + len(peaks)):
						if main.config["blendedlineswindow_autopositionpeaks"]:
							exp_mean = exp_ys.mean()

						yptp = 4*(np.amax(exp_ys)-np.amin(exp_ys))
						w0 = xrange[1] - xrange[0]
						wmax = main.config["blendedlineswindow_maxfwhm"] or w0
						amp_min, amp_max = -3*yptp, 3*yptp
						if amplitude_direction < 0:
							amp_max = 0
						if amplitude_direction > 0:
							amp_min = 0
						
						for peak in peaks:
							x, y, x_rel = peak

							if not xrange[0] < x < xrange[1]:
								x = self.center + x_rel
								if not xrange[0] < x < xrange[1]:
									x = sum(xrange)/2
								y = yptp * np.sign(amplitude_direction)
							elif not amp_min < y < amp_max:
								y = yptp * np.sign(amplitude_direction)

							xs.append((x, *xrange))
							ys.append((y, amp_min, amp_max))
							ws.append((wmax/4, 0, wmax))

						p0 = []
						bounds = [[], []]
						for x, y, w in zip(xs, ys, ws):
							tmp_p0 = [z[0] for z in (x, y, w, w)]
							tmp_b0 = [z[1] for z in (x, y, w, w)]
							tmp_b1 = [z[2] for z in (x, y, w, w)]

							if fixedwidth:
								slice_ = slice(2)
							elif profile == "Voigt":
								slice_ = slice(4)
							else:
								slice_ = slice(3)

							p0.extend(tmp_p0[slice_])
							bounds[0].extend(tmp_b0[slice_])
							bounds[1].extend(tmp_b1[slice_])

						if fixedwidth and len(peaks):
							ws = np.array(ws)
							w0, wl, wu = ws[:, 0].mean(), ws[:, 1].min(), ws[: 2].max()
							p0.extend([w0]*now)
							bounds[0].extend([wl]*now)
							bounds[1].extend([wu]*now)

						p0.extend([0]*polynomrank)
						bounds[0].extend([-np.inf]*polynomrank)
						bounds[1].extend([+np.inf]*polynomrank)

						# @Luis: Work on the bounds and p0s here
						try:
							popt, pcov = optimize.curve_fit(fitfunction, exp_xs, exp_ys, p0=p0, bounds=bounds)
						except Exception as E:
							popt, pcov = optimize.curve_fit(fitfunction, exp_xs, exp_ys, p0=p0, bounds=bounds)
						perr = np.sqrt(np.diag(pcov))
						res_xs = np.linspace(xrange[0], xrange[1], main.config["blendedlineswindow_xpoints"])
						res_ys = fitfunction(res_xs, *popt)
						res_exp_ys = fitfunction(exp_xs, *popt)
					else:
						popt = [0]*(noa*len(peaks)+polynomrank+2*fixedwidth)
						perr = [0]*(noa*len(peaks)+polynomrank+2*fixedwidth)
						res_xs = np.linspace(xrange[0], xrange[1], main.config["blendedlineswindow_xpoints"])
						res_ys = res_xs*0
						res_exp_ys = exp_xs*0

					earlyreturn(ownid, self.fit_thread_id)


					self.fit_line.set_data(res_xs, res_ys)
					self.fit_line.set_color(main.config["blendedlineswindow_color_total"])

					opt_param = []
					err_param = []

					for i in range(len(peaks)):
						tmp_params = list(popt[i*noa: (i+1)*noa])
						tmp_errors = list(perr[i*noa: (i+1)*noa])

						if fixedwidth:
							tmp_params.extend(popt[-(polynomrank+now):len(popt)-polynomrank])
							tmp_errors.extend(perr[-(polynomrank+now):len(popt)-polynomrank])

						tmp_ys = fitfunction_withoutbaseline(res_xs, *tmp_params)
						tmp_ys += exp_mean
						self.plot_parts.append(self.ax.plot(res_xs, tmp_ys, color=main.config["blendedlineswindow_color_total"], alpha=main.config["blendedlineswindow_transparency"])[0])

						opt_param.append( tmp_params )
						err_param.append( tmp_errors )


					self.plot_parts.append(self.ax.scatter([x[0] for x in opt_param], [x[1] + exp_mean for x in opt_param], color=main.config["blendedlineswindow_color_points"]))

					if polynomrank > 0 and main.config["blendedlineswindow_showbaseline"]:
						baseline_args = popt[-polynomrank:]
						self.plot_parts.append(self.ax.plot(res_xs, np.polyval(baseline_args, res_xs-self.center), color=main.config["blendedlineswindow_color_baseline"])[0])
					else:
						baseline_args = []

					rms_ys = np.sqrt(np.sum((exp_ys - res_exp_ys)**2) / len(exp_ys)) if len(exp_ys) else 0
					self.parent.params = opt_param, err_param, profile, derivative, noa, now, self.center, baseline_args, rms_ys

					earlyreturn(ownid, self.fit_thread_id)

					main.signalclass.blwfit.emit()
					main.signalclass.fitindicator.emit("Ready")
					super().update_plot()
				except CustomError as E:
					pass
				except Exception as E:
					main.signalclass.fitindicator.emit("<span style='color:#ff0000;'>Fit failed</span>")
					raise

			def update_plot(self):
				thread = threading.Thread(target=self.fit_peaks)
				with locks["currThread"]:
					thread.start()
					self.fit_thread_id = thread.ident
				return(thread)

		layout = QVBoxLayout()
		self.setLayout(layout)
		self.plotWidget = CustomPlotWidget(parent=self)
		layout.addWidget(self.plotWidget)

		tmplayout = QGridLayout()
		row_id = 0
		tmplayout.addWidget(QQ(QLabel, text="Function: "), row_id, 0)
		tmplayout.addWidget(QQ(QComboBox, "blendedlineswindow_lineshape", items=("Gauss", "Lorentz", "Voigt"), minWidth=120), row_id, 1)

		tmplayout.addWidget(QQ(QLabel, text="Transparency: "), row_id, 2)
		tmplayout.addWidget(QQ(QDoubleSpinBox, "blendedlineswindow_transparency", range=(0, 1), minWidth=120, singlestep=0.1), row_id, 3)

		row_id += 1
		tmplayout.addWidget(QQ(QLabel, text="Derivative: "), row_id, 0)
		tmplayout.addWidget(QQ(QSpinBox, "blendedlineswindow_derivative", range=(0, 2), minWidth=120), row_id, 1)

		row_id += 1
		tmplayout.addWidget(QQ(QLabel, text="Max FWHM: "), row_id, 0)
		tmplayout.addWidget(QQ(QDoubleSpinBox, "blendedlineswindow_maxfwhm", range=(0, None), minWidth=120), row_id, 1)

		tmplayout.addWidget(QQ(QCheckBox, "blendedlineswindow_fixedwidth", text="All Same Width"), row_id, 2, 1, 2)

		row_id += 1
		tmplayout.addWidget(QQ(QLabel, text="Baseline Rank: "), row_id, 0)
		tmplayout.addWidget(QQ(QSpinBox, "blendedlineswindow_polynom", range=(-1, None), minWidth=120), row_id, 1)

		tmplayout.addWidget(QQ(QCheckBox, "blendedlineswindow_showbaseline", text="Show Baseline"), row_id, 2, 1, 2)
		tmplayout.setColumnStretch(4, 1)

		layout.addLayout(tmplayout)

		self.label = QQ(QLabel, text="Ready", textFormat=Qt.TextFormat.RichText)
		main.signalclass.fitindicator.connect(self.label.setText)
		tmp_layout = QHBoxLayout()
		layout.addLayout(tmp_layout)
		row = (
		  QQ(QPushButton, text="Del All", change=lambda x: self.del_peak(-1)),
		  QQ(QPushButton, text="Update", change=lambda x: self.plotWidget.update_plot()),
		  QQ(QPushButton, text="Save", change=lambda x: self.save_values()),
		  self.label,
		)

		for widget in row:
			tmp_layout.addWidget(widget)
		tmp_layout.addStretch()

		self.table = QTableWidget()
		self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
		self.table.setMinimumHeight(50)
		layout.addWidget(self.table)

		self.cid = self.plotWidget.plot_canvas.mpl_connect("button_press_event", lambda event: self.add_peak(event.xdata, event.ydata))
		main.signalclass.blwfit.connect(self.fill_table)

	def add_peak(self, x, y):
		if not (x and y):
			return
		x = np.float64(x)
		y = np.float64(y)
		x_rel = x - self.plotWidget.center

		self.peaks.append((x, y, x_rel))
		self.peaks.sort(key=lambda x: x[0])
		self.plotWidget.update_plot()

	def del_peak(self, i=None):
		if i == -1:
			self.peaks = []
		elif i is not None and isinstance(i, (int, float, np.integer, np.float64)):
			if i in range(len(self.peaks)):
				del self.peaks[int(i)]
		else:
			if len(self.peaks) != 0:
				self.peaks.pop()
		self.plotWidget.update_plot()

	def fitfunction(self, x, fun, der, bpr, *ps, fixedwidth=False):
		noa = 2 + (1+(fun == "Voigt")) * (not fixedwidth)
		now = 2 if fun == "Voigt" else 1
		res_ys = []
		param_peaks = ps

		# Baseline Polynom
		if bpr > 0:
			param_baseline = ps[-bpr:]
			param_peaks    = ps[:-bpr]

			res_ys.append(np.polyval(param_baseline, x-self.plotWidget.center))

		# Fixed Width
		if fixedwidth:
			widths         = param_peaks[-now:]
			param_peaks    = param_peaks[:-now]

		# Peaks
		for i in range(len(param_peaks)//noa):
			tmp_params = list(param_peaks[i*noa: (i+1)*noa])
			if fixedwidth:
				tmp_params.extend(widths)
			res_ys.append(lineshape(fun, der, x, *tmp_params))

		return(np.sum(res_ys, axis=0))

	def fill_table(self):
		opt_param, err_param, function, derivative, noa, now, self.center, baseline_args, rms_ys = self.params
		fit_values = {
			"function": function,
			"derivative": derivative,
			"center": self.center,
			"baseline": list(baseline_args),
			"peaks": [],
			"RMS": rms_ys,
			"datetime": time.strftime("%d.%m.%Y %H:%M:%S", time.localtime()),
		}
		
		opt_param.sort(key=lambda x: x[0])
		table = self.table
		table.setRowCount(0)
		table.setColumnCount(7)
		table.setHorizontalHeaderLabels(["Action", "Frequency", "Amplitude", "FWHM Gauss", "FWHM Lorentz", "Other QNs", "Delete"])
		for params, params_error in zip(opt_param, err_param):
			x, y, wg, wl = params[0], params[1], params[2] if function != "Lorentz" else  0, params[1+now] if function != "Gauss" else  0
			x_error, y_error, wg_error, wl_error = params_error[0], params_error[1], params_error[2] if function != "Lorentz" else  0, params_error[1+now] if function != "Gauss" else  0
			
			currRowCount = table.rowCount()
			table.insertRow(currRowCount)
			table.setCellWidget(currRowCount, 0, QQ(QPushButton, text="Assign", change=lambda x, xpos=x, error=x_error: self.pre_assign(xpos, error)))
			table.setItem(currRowCount, 1, QTableWidgetItem(f'{x:{main.config["flag_xformatfloat"]}}'))
			table.setItem(currRowCount, 2, QTableWidgetItem(f'{y:{main.config["flag_xformatfloat"]}}'))
			table.setItem(currRowCount, 3, QTableWidgetItem(f'{wg:{main.config["flag_xformatfloat"]}}'))
			table.setItem(currRowCount, 4, QTableWidgetItem(f'{wl:{main.config["flag_xformatfloat"]}}'))
			table.setCellWidget(currRowCount, 5, QQ(QPushButton, text="Assign other QNs", change=lambda x, xpos=x, error=x_error: self.pre_assign(xpos, error, oqns=True)))
			table.setCellWidget(currRowCount, 6, QQ(QPushButton, text="Delete", change=lambda x, ind=currRowCount: self.del_peak(i=ind)))
			fit_values["peaks"].append({"values": (x, y, wg, wl), "errors": (x_error, y_error, wg_error, wl_error)})
		self.fit_values = fit_values
		table.resizeColumnsToContents()

	def save_values(self):
		if self.fit_values:
			filename = llwpfile(".fit")
			
			with open(filename, "a+") as file:
				file.seek(0)
				previous_fits = file.read()
				if previous_fits.strip():
					all_fits = json.loads(previous_fits)
					all_fits.append(self.fit_values)
				else:
					all_fits = [self.fit_values]
			
				file.truncate(0)
				json.dump(all_fits, file, indent=2)
			main.notification(f"Saved the fit to the file {filename}.")
		else:
			main.notification(f"No fit values to be saved.")

	def pre_assign(self, x, error, oqns=False):
		index = self.plotWidget.i

		dict_ = {
			"x": x,
			"error": error,
			"xpre": 0,
		}

		for i in range(6):
			dict_[f"qnu{i+1}"] = pyckett.SENTINEL
			dict_[f"qnl{i+1}"] = pyckett.SENTINEL

		if oqns:
			visible_cat_data = main.get_visible_data("cat")
			dialog = QNsDialog(x, visible_cat_data)
			dialog.exec()

			if dialog.result() == 1:
				dict_.update(dialog.save())
			else:
				return
		else:
			reference = main.plotwidget.get_series_reference(index[1])
			if reference["method"] == "Transition":
				qns = main.plotwidget.get_qns(reference["transition"])[index[0]]
				for i, qnu, qnl in zip(range(6), qns[0], qns[1]):
					dict_[f"qnu{i+1}"] = qnu
					dict_[f"qnl{i+1}"] = qnl
				dict_["xpre"] = main.plotwidget.cache_positions[index]

		if main.plotwidget.check_blends(index, dict_):
			return
		else:
			main.plotwidget.assign(index, dict_)

	def activateWindow(self):
		self.plotWidget.fig.canvas.mpl_disconnect(self.cid)
		self.cid = self.plotWidget.plot_canvas.mpl_connect("button_press_event", lambda event: self.add_peak(event.xdata, event.ydata))

		self.plotWidget.from_current_plot()
		super().activateWindow()

	@synchronized_d(locks["windows"])
	def closeEvent(self, *args, **kwargs):
		self.plotWidget.fig.canvas.mpl_disconnect(self.cid)
		return super().closeEvent(*args, **kwargs)

class SeriesfinderWindow(EQWidget):
	def __init__(self, id, parent=None):
		super().__init__(id, parent)
		self.setWindowTitle("Series Finder")

		layout = QVBoxLayout()

		self.messageLabel = QQ(QLabel, wordwrap=True, hidden=True)

		self.outputTable = QTableWidget()
		self.outputTable.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)

		vertLayout = QHBoxLayout()
		leftLayout = QGridLayout()
		rightLayout = QVBoxLayout()
		layout.addWidget(QQ(QLabel, wordwrap=True, text="Series Finder allows to find the strongest (unassigned) predicted transitions."))

		rightLayout.addWidget(QQ(QLabel, text="Allowed Transitions"))
		rightLayout.addWidget(QQ(QCheckBox, "seriesfinderwindow_atype", text="a-type"))
		rightLayout.addWidget(QQ(QCheckBox, "seriesfinderwindow_btype", text="b-type"))
		rightLayout.addWidget(QQ(QCheckBox, "seriesfinderwindow_ctype", text="c-type"))
		rightLayout.addStretch(1)

		leftLayout.addWidget(QQ(QLabel, text="Start Frequency: "), 1, 0)
		leftLayout.addWidget(QQ(QLineEdit, "seriesfinderwindow_start"), 1, 1)
		leftLayout.addWidget(QQ(QLabel, text="Stop Frequency: "), 2, 0, 1, 1)
		leftLayout.addWidget(QQ(QLineEdit, "seriesfinderwindow_stop"), 2, 1, 1, 1)
		leftLayout.addWidget(QQ(QLabel, text="Number of Results: "), 3, 0, 1, 1)
		leftLayout.addWidget(QQ(QSpinBox, "seriesfinderwindow_results", range=(0, None)), 3, 1, 1, 1)
		leftLayout.addWidget(QQ(QLabel, text="Additional Condition: "), 4, 0, 1, 1)
		leftLayout.addWidget(QQ(QLineEdit, "seriesfinderwindow_condition", placeholder="Hover for tooltip", tooltip="Use qnu1, ..., qnu6, qnl1, ..., qnl6 for the quantum numbers. Additionally, x, error, y, degfreed, elower, usd, tag, qnfmt, and filename are allowed."), 4, 1, 1, 1)
		leftLayout.addWidget(QQ(QCheckBox, "seriesfinderwindow_onlyunassigned", text = "Only unassigned Lines"), 5, 1, 1, 1)
		leftLayout.addWidget(QQ(QPushButton, text="Run", change=lambda x: self.run()), 6, 1, 1, 1)

		leftLayout.setColumnStretch(1, 1)

		vertLayout.addLayout(leftLayout)
		vertLayout.addLayout(rightLayout)
		layout.addLayout(vertLayout)
		layout.addWidget(self.messageLabel)
		layout.addWidget(self.outputTable, 1)
		self.setLayout(layout)

	@synchronized_d(locks["cat_df"])
	@synchronized_d(locks["lin_df"])
	def run(self):
		nor = main.config["seriesfinderwindow_results"]

		condition = []
		tmp_min = main.config["seriesfinderwindow_start"]
		tmp_max = main.config["seriesfinderwindow_stop"]
		addCondition = main.config["seriesfinderwindow_condition"]

		if addCondition.strip():
			condition.append(addCondition)
		if tmp_min:
			condition.append(f"{tmp_min} <= x")
		if tmp_max:
			condition.append(f"x <= {tmp_max}")

		tmp_condition = []
		if main.config["seriesfinderwindow_atype"]:
			tmp_condition.append(f"(abs(qnu2-qnl2) % 2 == 0 and abs(qnu3-qnl3) % 2 == 1)")
		if main.config["seriesfinderwindow_btype"]:
			tmp_condition.append(f"(abs(qnu2-qnl2) % 2 == 1 and abs(qnu3-qnl3) % 2 == 1)")
		if main.config["seriesfinderwindow_ctype"]:
			tmp_condition.append(f"(abs(qnu2-qnl2) % 2 == 1 and abs(qnu3-qnl3) % 2 == 0)")
		if tmp_condition:
			condition.append(" or ".join(tmp_condition))
		condition = " and ".join([f"({x})" for x  in condition])

		tmp_cat_df = main.get_visible_data("cat")

		if condition:
			try:
				tmp_cat_df.query(condition, inplace=True)
			except Exception as E:
				main.notification(f"<span style='color:#eda711;'>WARNING</span>: There is a syntax error in your condition: {str(E)}")
				return

		self.noq = noq = main.config["series_qns"]

		qns_visible = [f"qn{ul}{n+1}" for ul in ("u", "l") for n in range(noq)]
		qns_invisible = [f"qn{ul}{n+1}" for ul in ("u", "l") for n in range(noq, 6)]

		if main.config["seriesfinderwindow_onlyunassigned"]:
			tmp_lin_df = main.get_visible_data("lin")
			tmp_lin_df["DROP"] = True
			tmp_lin_df.drop(columns=["x"]+qns_invisible, inplace=True)

			tmp_cat_df = pd.merge(tmp_cat_df, tmp_lin_df, how="outer", on=qns_visible)
			tmp_cat_df = tmp_cat_df[tmp_cat_df.DROP != True]
			unassigned = "without already assigned lines"
		else:
			unassigned = "with already assigned lines"

		tmp_cat_df["y"] = np.log10(tmp_cat_df["y"])
		tmp_cat_df = tmp_cat_df.nlargest(nor, "y")

		if tmp_min and tmp_max:
			xrange = f"in the range from {tmp_min} to {tmp_max}"
		else:
			xrange = "in the total range"
		message = f"The {nor} most intense predicted transitions {unassigned} {xrange} are shown below."
		self.messageLabel.setText(message)
		self.messageLabel.setHidden(False)

		table = self.outputTable
		headers = ["Start", "Log. Intensity", "Frequency"] + qns_visible

		table.setRowCount(0)
		table.setColumnCount(len(headers))
		table.setHorizontalHeaderLabels(headers)

		for index, row in tmp_cat_df.iterrows():
			currRowCount = table.rowCount()
			table.insertRow(currRowCount)
			table.setCellWidget(currRowCount,0, QQ(QPushButton, text="Start", change=lambda x, crow=row: self.startHere(crow)))
			table.setItem(currRowCount, 1, QTableWidgetItem(f'{row["y"]:{main.config["flag_xformatfloat"]}}'))
			table.setItem(currRowCount, 2, QTableWidgetItem(f'{row["x"]:{main.config["flag_xformatfloat"]}}'))

			for i, column in enumerate(qns_visible):
				tmp = row[column]
				if tmp == pyckett.SENTINEL:
					tmp = ""
				else:
					tmp = f'{tmp:g}'
				table.setItem(currRowCount, i+3, QTableWidgetItem(tmp))

		self.outputTable.resizeColumnsToContents()

	def startHere(self, row):
		reference = main.config["series_references"][main.config["series_currenttab"]]
		reference["method"] = "Transition"
		reference["transition"]["qnus"][:self.noq] = [int(row[f"qnu{i+1}"]) for i in range(self.noq)]
		reference["transition"]["qnls"][:self.noq] = [int(row[f"qnl{i+1}"]) for i in range(self.noq)]
		main.config["series_qns"] = self.noq

		main.plotwidget.set_data()
		main.plotwidget.activateWindow()

class BlendWindow(EQWidget):
	def __init__(self, id, entries, dict_, index, parent=None):
		super().__init__(id, parent)
		self.setWindowTitle("Assign Blend")

		self.init_gui()
		self.update_gui(entries, dict_, index)

	def update_gui(self, entries, dict_, index):
		tmp_cat = main.get_visible_data("cat", (dict_["xpre"], dict_["xpre"]))
		qn_labels = [f"qn{ul}{i+1}" for ul in "ul" for i in range(6)]
		query = " and ".join([f"( {label} == {dict_[label]} )" for label in qn_labels if label in dict_])
		tmp_cat = tmp_cat.query(query)
		y_at_xpre = tmp_cat["y"].values[0]

		if main.config["series_blendminrelratio"] > 0:
			min_y_value = y_at_xpre * main.config["series_blendminrelratio"]
			entries = entries.query("y >= @min_y_value")

		self.noq = main.config["series_qns"]
		self.entries = entries.reset_index(drop=True)
		self.xmiddle = dict_["x"]
		self.index = index
		self.error = dict_["error"]

		self.label.setText(f"Assigning Blend for position {self.xmiddle}.")

		self.entries["dist"] = self.entries["x"] - self.xmiddle
		self.entries["absdist"] = abs(self.entries["dist"])
		self.entries.sort_values(by=["absdist"], inplace=True)
		self.entries["y"] = np.log(self.entries["y"])/np.log(10)

		noq = main.config["series_qns"]
		self.checkboxes = []
		table = self.table
		table.setRowCount(0)

		for i, entry in self.entries.iterrows():
			currRowCount = table.rowCount()
			table.insertRow(currRowCount)

			for j, val in enumerate(self.cols):
				if val not in entry:
					text = ""
				elif val == "filename":
					text = entry[val]
				else:
					text = f"{entry[val]:{main.config['flag_xformatfloat']}}".rstrip("0").rstrip(".")
				table.setItem(currRowCount, j+1, QTableWidgetItem(text))
			checkbox = QQ(QCheckBox, value=True)
			self.checkboxes.append(checkbox)
			table.setCellWidget(currRowCount, 0, checkbox)

		for i in range(6):
			table.setColumnHidden(i+ 4, i>=noq)
			table.setColumnHidden(i+10, i>=noq)

		table.resizeColumnsToContents()

	def init_gui(self):
		layout = QVBoxLayout()
		self.setLayout(layout)

		self.label = QQ(QLabel)
		layout.addWidget(self.label)

		self.table = QTableWidget()
		self.cols = ["x", "y", "dist"] + [f"qn{ul}{i+1}" for ul in ("u", "l") for i in range(6)] + ["filename"]
		self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
		self.table.setRowCount(0)
		self.table.setColumnCount(len(self.cols)+1)
		self.table.setHorizontalHeaderLabels(["Y/N", "x", "log. y", "Dist"] +  [f"{ul}{i+1}" for ul in ("U", "L") for i in range(6)] + ["Filename"])
		layout.addWidget(self.table)

		buttons_layout = QHBoxLayout()

		buttons_layout.addStretch()
		buttons_layout.addWidget(QQ(QPushButton, text="Assign All", shortcut="Ctrl+Return", change=lambda x: self.blendassign()))
		buttons_layout.addWidget(QQ(QPushButton, text="Assign Marked", shortcut="Return", change=lambda x: self.blendassign(False)))
		buttons_layout.addWidget(QQ(QPushButton, text="Select All", shortcut="Ctrl+A", change=lambda x: self.set_all(True)))
		buttons_layout.addWidget(QQ(QPushButton, text="Unselect All", shortcut="Ctrl+U", change=lambda x: self.set_all(False)))
		buttons_layout.addWidget(QQ(QPushButton, text="Cancel", change=lambda x: self.close()))
		buttons_layout.addStretch()

		layout.addLayout(buttons_layout)

	def set_all(self, state):
		for cb in self.checkboxes:
			cb.setChecked(state)

	def blendassign(self, all=True):
		if all == True:
			checked = True
		else:
			checked = [x.isChecked() for x in self.checkboxes]
		self.entries["checked"] = checked

		new_df = main.new_df

		selected_rows = self.entries.query("checked == True")
		max_predicted = selected_rows["y"].max()

		for i, row in selected_rows.iterrows():
			tmp_dict = {
				"x": 		self.xmiddle,
				"weight":	10 ** (row["y"]-max_predicted),
				"error":	self.error,
				"comment":	main.config["fit_comment"],
				"filename":	"__lin__",
			}

			for i in range(6):
				tmp_dict[f"qnu{i+1}"] = pyckett.SENTINEL
				tmp_dict[f"qnl{i+1}"] = pyckett.SENTINEL

			for i in range(self.noq):
				tmp_dict[f"qnu{i+1}"] = row[f"qnu{i+1}"]
				tmp_dict[f"qnl{i+1}"] = row[f"qnl{i+1}"]
			main.plotwidget.assign(None, tmp_dict)

		main.plotwidget.set_data()
		self.close()

class PipeWindow(EQWidget):
	def __init__(self, id, parent=None):
		super().__init__(id, parent)
		self.setWindowTitle("Pipe")

		self.tabs = QTabWidget()
		self.tabs.setTabsClosable(True)
		self.tabs.setMovable(True)
		self.tabs.setDocumentMode(True)

		if not main.config["pipe_commands"]:
			main.config["pipe_commands"] = [["Initial", "", True, True, True]]

		for values in main.config["pipe_commands"]:
			self.add_tab(*values)

		self.tabs.tabCloseRequested.connect(self.close_tab)
		self.tabs.tabBarDoubleClicked.connect(self.renameoradd_tab)
		self.tabs.setCurrentIndex(main.config["pipe_current"])
		self.tabs.currentChanged.connect(lambda x: self.update_pipe_command())

		layout = QVBoxLayout()
		layout.addWidget(QQ(QLabel, wordwrap=True, text="Here you can set up a command to run fitting and prediction programs to update your files."))
		layout.addWidget(self.tabs)
		layout.addWidget(QQ(QPushButton, text="Run", change=lambda x: main.run_pipe()))
		self.setLayout(layout)

	def add_tab(self, title, command=None, exp_checked=True, cat_checked=True, lin_checked=True):
		tmp = QWidget()
		layout = QVBoxLayout()

		if not command:
			command = ""

		layout.addWidget(QQ(QPlainTextEdit, value=command, change=lambda: self.update_pipe_command(), placeholder="Write your command line command here"))
		layout.addWidget(QQ(QCheckBox, text="Reread Exp Files after command finished", value=exp_checked, change=lambda x: self.update_pipe_command()))
		layout.addWidget(QQ(QCheckBox, text="Reread Cat Files after command finished", value=cat_checked, change=lambda x: self.update_pipe_command()))
		layout.addWidget(QQ(QCheckBox, text="Reread Lin Files after command finished", value=lin_checked, change=lambda x: self.update_pipe_command()))

		tmp.setLayout(layout)
		self.tabs.addTab(tmp, title)

	@synchronized_d(locks["pipe"])
	def close_tab(self, index):
		tab = self.tabs.widget(index)
		tab.deleteLater()
		self.tabs.removeTab(index)
		if len(main.config["pipe_commands"]) > index:
			main.config["pipe_commands"].pop(index)
		if self.tabs.count() == 0:
			self.add_tab("New Tab")
			main.config["pipe_commands"].append(["New Tab", "", True, True, True])

	@synchronized_d(locks["pipe"])
	def renameoradd_tab(self, index):
		if index == -1:
			self.add_tab("New Tab")
			main.config["pipe_commands"].append(["New Tab", "", True, True, True])
		elif self.tabs.widget(index) != 0:
			text, ok = QInputDialog().getText(self, "Tab Name","Enter the Tabs Name:")
			if ok and text:
				self.tabs.setTabText(index, text)
				main.config["pipe_commands"][index][0] = text

	@synchronized_d(locks["pipe"])
	def update_pipe_command(self):
		result = []
		for i in range(self.tabs.count()):
			tab = self.tabs.widget(i)
			title = self.tabs.tabText(i)
			command = tab.findChildren(QPlainTextEdit)[0].toPlainText()
			if command.strip() == "":
				command = None
			rereadCheckboxs = tab.findChildren(QCheckBox)
			reread_files = [x.isChecked() for x in rereadCheckboxs]

			tmp = [title, command, *reread_files]
			result.append(tmp)

		main.config["pipe_commands"] = result
		main.config["pipe_current"] = self.tabs.currentIndex()

class PeakfinderWindow(EQWidget):
	@synchronized_d(locks["peaks"])
	def __init__(self, id, parent=None):
		super().__init__(id, parent)
		self.setWindowTitle("Peakfinder")

		self.peaks = []

		class CustomPlotWidget(ProtPlot):
			def gui(self):
				super().gui()
				self.peaks_line = self.ax.scatter([], [], color=main.config["peakfinderwindow_peakcolor"], marker="*")

			def update_plot(self):
				xmin, xmax = self.center-self.width/2, self.center+self.width/2
				tuples = list(filter(lambda x: xmin < x[0] < xmax, self.parent.peaks))
				tuples = tuples if len(tuples)!=0 else [[None,None]]
				self.peaks_line.set_offsets(tuples)
				self.peaks_line.set_color(main.config["peakfinderwindow_peakcolor"])
				super().update_plot()

		layout = QVBoxLayout()
		self.setLayout(layout)

		self.plotWidget = CustomPlotWidget(parent=self)
		layout.addWidget(self.plotWidget)

		tmp_layout = QHBoxLayout()
		layout.addLayout(tmp_layout)

		self.run_button = QQ(QPushButton, text="Run Peakfinder", change=lambda: self.find_peaks())
		tmp_layout.addWidget(self.run_button)
		tmp_layout.addWidget(QQ(QPushButton, text="Export Peaks", change=lambda: self.export_peaks()))

		uncert_input = QQ(QDoubleSpinBox, "peakfinderwindow_width", range=(0, None), enabled=main.config["peakfinderwindow_onlyunassigned"], minWidth=80)
		tmp_layout.addWidget(QQ(QCheckBox, "peakfinderwindow_onlyunassigned", text="Only unassigned Lines", change=uncert_input.setEnabled))
		tmp_layout.addWidget(QQ(QLabel, text="Uncertainty unassigned lines: ", tooltip="Max distance between peak and closest predicted line to count as assigned"))
		tmp_layout.addWidget(uncert_input)

		tmp_layout.addStretch(1)

		tmp_layout = QGridLayout()
		layout.addLayout(tmp_layout)
		tmp_layout.addWidget(QQ(QLabel, text="Configure Peakfinder:"), 0, 0)
		tmp_layout.addWidget(QQ(QLabel, text="Min:"), 0, 1)
		tmp_layout.addWidget(QQ(QLabel, text="Max:"), 0, 2)

		self.row_dict = {}
		for i, key in enumerate(("frequency", "height", "threshold", "distance", "prominence", "width")):
			label      = QQ(QLabel, text=key.capitalize())
			input_min  = QQ(QDoubleSpinBox, range=(0, None))
			input_max  = QQ(QDoubleSpinBox, range=(0, None))

			tmp_layout.addWidget(label, i+1, 0)
			tmp_layout.addWidget(input_min, i+1, 1)
			if key == "distance":
				items = (("Min", "Off"))
				input_min.setRange(1, None)
			else:
				items = (("Min", "Min & Max", "Off"))
				tmp_layout.addWidget(input_max, i+1, 2)
			input_type = QQ(QComboBox, items=items, change=lambda x, key=key: self.change_type(key))
			tmp_layout.addWidget(input_type, i+1, 3)

			self.row_dict[key] = (input_min, input_max, input_type, label)

			if key in main.config["peakfinderwindow_kwargs"]:
				value = main.config["peakfinderwindow_kwargs"][key]
				try:
					if type(value) in [tuple, list]:
						input_min.setValue(value[0])
						input_max.setValue(value[1])
						input_type.setCurrentIndex(1)
					else:
						input_min.setValue(value)
						input_type.setCurrentIndex(0)
				except:
					pass
			else:
				input_type.setCurrentIndex(len(items)-1)
			self.change_type(key)

		self.infolabel = QQ(QLabel, wordwrap=True, text="Press 'Run Peakfinder' to find peaks.")
		layout.addWidget(self.infolabel)

		self.table = QTableWidget()
		self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
		self.table.setHidden(True)
		self.table.doubleClicked.connect(self.go_to)

		layout.addWidget(self.table)

		main.signalclass.peakfinderstart.connect(lambda: self.run_button.setEnabled(False))
		main.signalclass.peakfinderstart.connect(lambda: self.infolabel.setText("Finding peaks ..."))
		main.signalclass.peakfinderend.connect(self.after_find_peaks)

	def change_type(self, key):
		input_min, input_max, input_type, label = self.row_dict[key]
		type = input_type.currentText()

		if type == "Min":
			enabled = (True, False)
		elif type == "Min & Max":
			enabled = (True, True)
		else:
			enabled = (False, False)

		input_min.setEnabled(enabled[0])
		input_max.setEnabled(enabled[1])
		label.setEnabled(enabled[0] or enabled[1])

	def get_kwargs(self):
		kwargs = {}
		for key, value in self.row_dict.items():
			input_min, input_max, input_type, label = value
			type = input_type.currentText()

			if type == "Min":
				kwargs[key] = input_min.value()
			elif type == "Min & Max":
				kwargs[key] = (input_min.value(), input_max.value())
		main.config["peakfinderwindow_kwargs"] = kwargs
		return(kwargs.copy())

	@threading_d
	@synchronized_d(locks["peaks"])
	def find_peaks(self):
		main.signalclass.peakfinderstart.emit()
		kwargs = self.get_kwargs()

		with locks["exp_df"]:
			if "frequency" in kwargs:
				val = kwargs["frequency"]
				if type(val) in [tuple, list]:
					exp_df = main.get_visible_data("exp", xrange=val)
				else:
					exp_df = main.get_visible_data("exp").query(f"{val} < x")
				del kwargs["frequency"]
			else:
				exp_df = main.get_visible_data("exp")
			xs, ys = exp_df["x"].to_numpy(), exp_df["y"].to_numpy()

		peaks, props = signal.find_peaks(ys, **kwargs)
		self.peaks = np.array((xs[peaks], ys[peaks])).T
		self.plotWidget.update_plot()

		if main.config["peakfinderwindow_onlyunassigned"]:
			assigned_xs = main.get_visible_data("lin")["x"]
			uncertainty = main.config["peakfinderwindow_width"]

			peaks_xs = self.peaks[:, 0]
			peaks_assigned = np.zeros(self.peaks.shape[0])

			for x in assigned_xs:
				peaks_assigned += (abs(peaks_xs - x) < uncertainty)
			self.peaks = self.peaks[peaks_assigned == 0]

		self.peaks = self.peaks[self.peaks[:, 1].argsort()[::-1]]
		main.signalclass.peakfinderend.emit()

	@synchronized_d(locks["peaks"])
	def after_find_peaks(self):
		self.table.setRowCount(0)
		self.table.setColumnCount(2)
		self.table.setHorizontalHeaderLabels(["x", "y"])
		for x, y in self.peaks[:main.config["peakfinderwindow_maxentries"]]:
			currRowCount = self.table.rowCount()
			self.table.insertRow(currRowCount)
			self.table.setItem(currRowCount, 0, QTableWidgetItem(f'{x:{main.config["flag_xformatfloat"]}}'))
			self.table.setItem(currRowCount, 1, QTableWidgetItem(f'{y:{main.config["flag_xformatfloat"]}}'))
		self.table.resizeColumnsToContents()
		self.table.setHidden(False)

		self.run_button.setEnabled(True)

		if len(self.peaks) > main.config["peakfinderwindow_maxentries"]:
			self.infolabel.setText(f"Found {len(self.peaks)} peaks. Only the highest {main.config['peakfinderwindow_maxentries']} entries are shown in the table for better performance. You can see all peaks by exporting or increasing the maximum number of displayed peaks in the configuration.")
		else:
			self.infolabel.setText(f"Found {len(self.peaks)} peaks.")

	@synchronized_d(locks["peaks"])
	def export_peaks(self):
		fname = QFileDialog.getSaveFileName(None, 'Choose file to save peaks to',"","CSV Files (*.csv);;All Files (*)")[0]
		if fname:
			np.savetxt(fname, self.peaks, delimiter="\t")

	@synchronized_d(locks["peaks"])
	def go_to(self, event):
		for idx in self.table.selectionModel().selectedIndexes():
			row_number = idx.row()

		if type(row_number) == int and len(self.peaks) > row_number:
			self.plotWidget.center = self.peaks[row_number, 0]
			self.plotWidget.update_plot()

	def activateWindow(self):
		self.plotWidget.from_current_plot()
		super().activateWindow()

class SeriesFitWindow(EQWidget):
	def __init__(self, id, parent=None):
		super().__init__(id, parent)
		self.setWindowTitle("Series Fit")

		self.pred_function = None
		self.pred_params = {}
		self.pred_freqs = []

		self.catalogTable = QTableView()
		self.catalogModel = CustomTableModel(main.ser_df, list(main.ser_df.columns))
		self.catalogTable.setModel(self.catalogModel)
		self.catalogTable.resizeColumnsToContents()

		class CustomSeriesSelector(SeriesSelector):
			def changed(self):
				super().changed()
				main.config["seriesfitwindow_series"] = self.state

		self.series_selector = CustomSeriesSelector(self, main.config["seriesfitwindow_series"])
		self.output_area = QQ(QPlainTextEdit, readonly=True, minHeight=50, placeholder="This is the output field and is read-only.")

		layout = QVBoxLayout()
		buttonsBox = QHBoxLayout()
		buttonsBox.addWidget(QQ(QToolButton, text="X", change=lambda x: self.catalogTableDelete()))
		buttonsBox.addWidget(QQ(QToolButton, text="+", change=lambda x: addemptyrow_inplace(main.ser_df, self.catalogModel)))
		buttonsBox.addWidget(QQ(QToolButton, text="Copy Lines From Main Window", change=lambda x: self.copyMWLines()))
		buttonsBox.addWidget(QQ(QCheckBox, "seriesfitwindow_greedy", text="Add Newly Assigned Lines here"))
		buttonsBox.addStretch(1)
		buttonsBox.addWidget(QQ(QToolButton, text="Tight", change=lambda x: self.catalogTable.resizeColumnsToContents()))
		layout.addLayout(buttonsBox)
		layout.addWidget(self.catalogTable)
		layout.addWidget(self.series_selector)
		layout.addWidget(QQ(QPlainTextEdit, "seriesfitwindow_function", minHeight=50, placeholder="Write formula for the fit here. The formula can include the six upper quantum numbers (as qnu1, ..., qnu6) and six lower quantum numbers (as qnl1, ..., qnl6). E.g. try 2*B*qnu1 for a molecule with equidistant peaks of about 8000 units."))

		hbox = QHBoxLayout()
		hbox.addWidget(QQ(QPushButton, text="Fit", change=lambda x: self.fit()))
		hbox.addWidget(QQ(QPushButton, text="Show", change=lambda x: self.show_pred_freqs()))
		hbox.addStretch(1)

		layout.addLayout(hbox)
		layout.addWidget(self.output_area)
		self.setLayout(layout)

		main.signalclass.assignment.connect(self.scroll_bottom)

	def scroll_bottom(self):
		self.catalogTable.selectRow(len(self.catalogModel.data)-1)
		self.catalogModel.update()
		self.catalogTable.scrollToBottom()

	@synchronized_d(locks["ser_df"])
	def catalogTableDelete(self):
		selected = [index.row() for index in self.catalogTable.selectionModel().selectedRows()]
		for index in sorted(selected, reverse = True):
			main.ser_df.drop(index, inplace =True)
			main.ser_df.reset_index(inplace = True, drop = True)
		self.catalogModel.update()

	def copyMWLines(self):
		main.ser_df.reset_index(drop=True, inplace=True)
		main.new_df.reset_index(drop=True, inplace=True)
		length = len(main.ser_df)
		tmp_df = main.new_df[main.ser_df.columns]

		for i in range(len(main.new_df)):
			tmp_dict = tmp_df.loc[i]
			tmp_dict = {key: main.ser_dtypes[key](value) for key, value in tmp_dict.items()}
			main.ser_df.loc[length+i] = tmp_dict
		self.catalogModel.update()

	def fit(self):
		state = main.config["seriesfitwindow_series"]
		self.noq = noq = main.config["series_qns"]
		cols = [f"qn{ul}{i+1}" for ul in ("u", "l") for i in range(noq)]

		incr_values = state["qndiffs"] if state["increase_freely"] else state["qnincrs"]
		qnus = state["qnus"]
		qnls = state["qnls"]

		conditions = []
		conditions_incr = []
		for i, qnu, qnl, incr in zip(range(noq), qnus, qnls, incr_values):
			if incr:
				conditions_incr.append(f"(qnu{i+1} - {qnu})*{incr}")
				conditions_incr.append(f"(qnl{i+1} - {qnl})*{incr}")
			else:
				conditions.append(f"(qnu{i+1} == {qnu})")
				conditions.append(f"(qnl{i+1} == {qnl})")

		if len(conditions_incr):
			conditions.append(" == ".join(conditions_incr))

		conditions = " and ".join(conditions)

		df = main.ser_df.query(conditions).copy()
		df.sort_values("x", inplace=True)

		if not len(df):
			self.writelog("No lines matching the transition.")
			return

		fitfunction = main.config["seriesfitwindow_function"]
		if not fitfunction:
			self.writelog("No fitfunction specified.")
			return

		fitfunction = compile(fitfunction, "<string>", "eval")
		fitparams   = set(fitfunction.co_names)-set(cols)

		if len(df) < len(fitparams):
			self.writelog("Too few assignments for too many parameters.")
			return
		elif len(fitparams) == 0:
			self.writelog("No free parameters in your expression.")
			return


		p0 = [1]*len(fitparams)
		self.fitparams = fitparams
		self.fitfunction = fitfunction
		try:
			popt, pcov = optimize.curve_fit(self.function, df[cols].to_numpy().transpose(), df["x"].to_numpy(), p0=p0)
		except Exception as E:
			popt, pcov = optimize.curve_fit(self.function, df[cols].to_numpy().transpose(), df["x"].to_numpy(), p0=p0)

		qnus = np.array(qnus[:noq])
		qnls = np.array(qnls[:noq])
		incrs = np.array(incr_values[:noq])

		self.pred_qns = np.array([np.concatenate((qnus+incrs*i, qnls+incrs*i)) for i in range(main.config["seriesfitwindow_maxprediction"])])
		self.pred_xs  = [self.function(qns, *popt) for qns in self.pred_qns]

		tmp = "\n".join([f"{name} : {value:{main.config['flag_xformatfloat']}}" for name, value in zip(self.fitparams, popt)])
		self.writelog(f"Succeeded, parameters were determined as \n{tmp}")

	def show_pred_freqs(self):
		original_shape = self.pred_qns.shape
		new_shape = (original_shape[0], 2, -1)
		qns = self.pred_qns.reshape(new_shape).tolist()
		reference = main.config["series_references"][main.config["series_currenttab"]]
		reference["method"] = "List"
		reference["list"] = {"qns": qns, "xs": self.pred_xs, "i0": 0}
		main.config["series_qns"] = self.noq

		main.mainwindow.referenceserieswindow.set_state(main.config["series_references"])
		main.plotwidget.set_data()
		main.plotwidget.activateWindow()

	def writelog(self, text):
		time_str = time.strftime("%H:%M", time.localtime())
		self.output_area.appendPlainText(f"{time_str}: {text}\n")
		sb = self.output_area.verticalScrollBar()
		sb.setValue(sb.maximum())

	def function(self, qns, *params):
		keys = [f"qn{ul}{i+1}" for ul in ("u", "l") for i in range(self.noq)]
		params_dict = {key: value for key, value in zip(keys, qns)}

		for name, param in zip(self.fitparams, params):
			params_dict[name] = param

		return(eval(self.fitfunction, params_dict))

class ResidualsWindow(EQWidget):
	def __init__(self, id, parent=None):
		super().__init__(id, parent)
		self.setWindowTitle("Residuals")

		self.fig = figure.Figure(dpi=main.config["plot_dpi"])
		self.plot_canvas = FigureCanvas(self.fig)

		self.ax = self.fig.subplots()
		self.ax.ticklabel_format(useOffset=False)

		self.points = self.ax.scatter([], [], color=[], marker=".")
		self.annot = self.ax.annotate("", xy=(0,0), xytext=(5,5), textcoords="offset points", color="black", ha="center", va="bottom", bbox=dict(boxstyle="round", fc="w"))
		self.annot.set_visible(False)
		self.fig.canvas.mpl_connect("motion_notify_event", self.on_hover)
		self.fig.canvas.mpl_connect("button_press_event", lambda x: self.on_hover(x, True))

		self.mpl_toolbar = NavigationToolbar2QT(self.plot_canvas, self)
		self.df = None

		layout = QVBoxLayout()
		self.setLayout(layout)
		layout.addWidget(self.plot_canvas, 6)
		hlayout = QHBoxLayout()
		hlayout.addWidget(QLabel("x-axis: "))
		hlayout.addWidget(QQ(QLineEdit, "residualswindow_xvariable", placeholder="Choose the x-variable, e.g. x_lin"))
		hlayout.addWidget(QLabel("y-axis: "))
		hlayout.addWidget(QQ(QLineEdit, "residualswindow_yvariable", placeholder="Choose the y-variable, e.g. x_lin-x_cat"))
		hlayout.addWidget(QQ(QCheckBox, "residualswindow_blends", text="Blends"))
		hlayout.addWidget(QQ(QCheckBox, "residualswindow_autoscale", text="Autoscale on Update"))
		layout.addLayout(hlayout)
		layout.addWidget(self.mpl_toolbar)
		layout.addWidget(QQ(QPlainTextEdit, "residualswindow_query", maxHeight=60, placeholder="Query text to filter shown lines. Use qnu1, ..., qnu6 and qnl1, ..., qnl6 for the quantum numbers. Other possible values are x_lin, x_cat, error_lin, error_cat, filename_lin, filename_cat, y, degfreed, elower, usd, tag, qnfmt, weight, and comment."))
		layout.addWidget(QQ(QPlainTextEdit, "residualswindow_colorinput", maxHeight=60, placeholder="Enter custom color and query to color specific lines differently. E.g. enter '#ff0000; qnu1 < 20' to color all transitions with the first upper quantum number below 20 red."))

		buttonslayout = QHBoxLayout()
		layout.addLayout(buttonslayout)
		buttonslayout.addStretch(1)
		self.update_button = QQ(QPushButton, text="Update", change=self.plot_residuals)
		buttonslayout.addWidget(self.update_button)
		buttonslayout.addWidget(QQ(QPushButton, text="Save", change=self.save_residuals))
		buttonslayout.addStretch(1)

	def get_residuals(self):
		lin_df = main.get_visible_data("lin")
		cat_df = main.get_visible_data("cat")

		noq = main.config["series_qns"]
		self.noq = noq
		qns_visible = [f"qn{ul}{n+1}" for ul in ("u", "l") for n in range(noq)]
		df = pd.merge(lin_df, cat_df, how="inner", on=qns_visible)
		df.rename(columns={"x_x": "x_lin", "x_y": "x_cat", "error_x": "error_lin", "error_y": "error_cat", "filename_x": "filename_lin", "filename_y": "filename_cat"}, inplace=True)

		query = main.config["residualswindow_query"].strip()
		if query:
			df.query(query, inplace=True)
		df.reset_index(drop=True, inplace=True)

		if main.config["residualswindow_blends"]:
			mask = df["weight"] != 0
			df_view = df[mask]
			tmp_dict = df_view.groupby(df_view.x_lin).apply(lambda x: np.average(x.x_cat, weights=x.weight)).to_dict()
			df["x_cat"] = df["x_lin"].map(lambda x: tmp_dict.get(x, x))

		df["obs_calc"] = df["x_lin"] - df["x_cat"]
		return(df)

	def plot_residuals(self):
		self.update_button.setDisabled(True)
		main.app.processEvents()
		try:
			message = []

			df = self.get_residuals()
			df["color"] = main.config["residualswindow_defaultcolor"]
			message.append(f"Found {len(df)} entries matching your query.")

			colorquerytext = main.config["residualswindow_colorinput"].split("\n")
			for row in colorquerytext:
				if row.strip():
					color, query = row.split(";")
					indices = df.query(query).index
					df.loc[indices, "color"] = color
					message.append(f"{len(indices)} lines are colored in <span style='color:{color};'>{color}</span>.")

			self.df = df

			yvariable = main.config["residualswindow_yvariable"].strip() or "obs_calc"
			xvariable = main.config["residualswindow_xvariable"].strip() or "x_lin"
			ys = df.eval(yvariable).to_numpy()
			xs = df.eval(xvariable).to_numpy()
			colors = df["color"].to_numpy()
			tuples = list(zip(xs,ys))
			tuples = tuples if len(tuples)!=0 else [[None,None]]
			self.points.set_offsets(tuples)
			self.points.set_color(colors)
			if len(xs) and main.config["residualswindow_autoscale"]:
				xmin, xmax = np.min(xs), np.max(xs)
				if xmin == xmax:
					xmin -= 1
					xmax += 1
				self.ax.set_xlim([xmin, xmax])
				ymin, ymax = np.min(ys), np.max(ys)
				if ymin == ymax:
					ymin -= 1
					ymax += 1
				y_range = [ymin, ymax]
				self.ax.set_ylim(y_range[0]-main.config["plot_ymargin"]*(y_range[1]-y_range[0]), y_range[1]+main.config["plot_ymargin"]*(y_range[1]-y_range[0]))
			self.fig.canvas.draw()
			main.notification("<br/>".join(message))
		except:
			main.notification("<span style='color:#eda711;'>WARNING</span>: There was an error in your Residuals window input")
			raise
		finally:
			self.update_button.setDisabled(False)

	def save_residuals(self):
		df = self.get_residuals()
		fname = QFileDialog.getSaveFileName(None, 'Choose file to save residuals to',"","CSV Files (*.csv);;All Files (*)")[0]
		if fname:
			df.to_csv(fname, index=None, sep='\t')

	def on_hover(self, event, click=False):
		if event.inaxes == self.ax and isinstance(self.df, pd.DataFrame):
			cont, ind = self.points.contains(event)
			noq = None
			if cont:
				self.annot.xy = self.points.get_offsets()[ind["ind"][0]]
				tmp_transitions = self.df.iloc[ind["ind"]]
				text = []
				if len(tmp_transitions):
					noq = self.noq
				for i, row in tmp_transitions.iterrows():
					qnus = [row[f"qnu{n+1}"] for n in range(noq)]
					qnls = [row[f"qnl{n+1}"] for n in range(noq)]
					text.append(f"{', '.join(map(str, qnus))} ← {', '.join(map(str, qnls))}")
				text = "\n".join(text)
				self.annot.set_text(text)
				self.annot.set_visible(True)
			else:
				self.annot.set_visible(False)
			self.fig.canvas.draw()

			if click and noq:
				# @Luis: can this be out of range
				reference = main.config["series_references"][main.config["series_currenttab"]]
				reference["method"] = "Transition"
				reference["transition"]["qnus"][:noq] = qnus
				reference["transition"]["qnls"][:noq] = qnls
				main.config["series_qns"] = noq

				main.plotwidget.set_data()
				main.plotwidget.activateWindow()

class EnergyLevelsWindow(EQWidget):
	def __init__(self, id, parent=None):
		super().__init__(id, parent)
		self.setWindowTitle("Energy Levels")
		self.setAcceptDrops(True)

		self.fig = figure.Figure(dpi=main.config["plot_dpi"])
		self.plot_canvas = FigureCanvas(self.fig)

		self.fname = None
		self.dataframe = None
		self.dataframe_filtered = None

		self.ax = self.fig.subplots()
		self.ax.ticklabel_format(useOffset=False)

		self.points = self.ax.scatter([], [], color=[], marker=".")
		self.annot = self.ax.annotate("", xy=(0,0), xytext=(5,5), textcoords="offset points", color="black", ha="center", va="bottom", bbox=dict(boxstyle="round", fc="w"))
		self.annot.set_visible(False)
		self.fig.canvas.mpl_connect("motion_notify_event", self.on_hover)

		self.mpl_toolbar = NavigationToolbar2QT(self.plot_canvas, self)

		layout = QVBoxLayout()
		self.setLayout(layout)
		layout.addWidget(self.plot_canvas, 6)
		layout.addWidget(self.mpl_toolbar)
		hlayout = QHBoxLayout()
		hlayout.addWidget(QQ(QPushButton, text="Open", change=lambda x: self.load_file()))
		self.file_label = QQ(QLabel, text="No File loaded")
		hlayout.addWidget(self.file_label)
		hlayout.addWidget(QLabel("x-axis: "))
		hlayout.addWidget(QQ(QLineEdit, "energylevelswindow_xvariable", placeholder="Choose the x-variable, e.g. qn1"))
		hlayout.addWidget(QLabel("y-axis: "))
		hlayout.addWidget(QQ(QLineEdit, "energylevelswindow_yvariable", placeholder="Choose the x-variable, e.g. egy"))
		hlayout.addWidget(QQ(QCheckBox, "energylevelswindow_autoscale", text="Autoscale on Update"))
		layout.addLayout(hlayout)
		layout.addWidget(QQ(QPlainTextEdit, "energylevelswindow_query", maxHeight=40, placeholder="Query text to filter shown levels. Use qn1, ..., qn6 for the quantum numbers. Other possible values are iblk, indx, egy, err, pmix, and we."))
		layout.addWidget(QQ(QPlainTextEdit, "energylevelswindow_colorinput", maxHeight=40, placeholder="Enter custom color and query to color specific lines differently. E.g. enter '#ff0000; qn1 < 20' to color all levels with the first quantum number below 20 red."))

		buttonslayout = QHBoxLayout()
		layout.addLayout(buttonslayout)
		buttonslayout.addStretch(1)
		self.update_button = QQ(QPushButton, text="Update", change=self.plot_energylevels)
		buttonslayout.addWidget(self.update_button)
		buttonslayout.addStretch(1)

	def load_file(self, fname=None):
		if fname is None:
			fname = QFileDialog.getOpenFileName(None, 'Choose Egy File to load',"")[0]
		if fname:
			self.dataframe = pyckett.egy_to_df(fname)
			self.dataframe_filtered = self.dataframe
			self.fname = fname
			self.file_label.setText(os.path.split(fname)[1])
			self.plot_energylevels()

	def plot_energylevels(self):
		self.update_button.setDisabled(True)
		main.app.processEvents()
		self.noq = main.config["series_qns"]
		try:
			if self.dataframe is None:
				return
			df = self.dataframe
			query = main.config["energylevelswindow_query"]
			if query:
				df = df.query(query).copy()

			df.loc[:, "color"] = main.config["energylevelswindow_defaultcolor"]
			colorquerytext = main.config["energylevelswindow_colorinput"].split("\n")
			for row in colorquerytext:
				if row.strip():
					color, query = row.split(";")
					df.loc[df.query(query).index, "color"] = color

			self.dataframe_filtered = df

			xvariable = main.config["energylevelswindow_xvariable"].strip() or "qn1"
			yvariable = main.config["energylevelswindow_yvariable"].strip() or "egy"
			xs = df.eval(xvariable).to_numpy()
			ys = df.eval(yvariable).to_numpy()
			colors = df["color"].to_numpy()
			tuples = list(zip(xs,ys))
			tuples = tuples if len(tuples)!=0 else [[None,None]]
			self.points.set_offsets(tuples)
			self.points.set_color(colors)
			if len(xs) and main.config["energylevelswindow_autoscale"]:
				xmin, xmax = np.min(xs), np.max(xs)
				if xmin == xmax:
					xmin -= 1
					xmax += 1
				self.ax.set_xlim([xmin, xmax])
				ymin, ymax = np.min(ys), np.max(ys)
				if ymin == ymax:
					ymin -= 1
					ymax += 1
				y_range = [ymin, ymax]
				self.ax.set_ylim(y_range[0]-main.config["plot_ymargin"]*(y_range[1]-y_range[0]), y_range[1]+main.config["plot_ymargin"]*(y_range[1]-y_range[0]))
			self.fig.canvas.draw()
		except:
			main.notification("There was an error in your Energy Levels window inputs")
			raise
		finally:
			self.update_button.setDisabled(False)
		

	def on_hover(self, event):
		if event.inaxes == self.ax and isinstance(self.dataframe_filtered, pd.DataFrame):
			cont, ind = self.points.contains(event)
			if cont:
				self.annot.xy = self.points.get_offsets()[ind["ind"][0]]
				tmp_levels = self.dataframe_filtered.iloc[ind["ind"]]
				text = []
				for i, row in tmp_levels.iterrows():
					text.append(",".join(str(int(row[f"qn{i+1}"])) for i in range(self.noq)))
				text = "\n".join(text)
				self.annot.set_text(text)
				self.annot.set_visible(True)
			else:
				self.annot.set_visible(False)
			self.fig.canvas.draw()

	def dragEnterEvent(self, event):
		mimeData = event.mimeData()
		if mimeData.hasUrls() and len(mimeData.urls()) == 1:
			event.accept()
		else:
			event.ignore()

	def dropEvent(self, event):
		url = event.mimeData().urls()[0]
		fname = url.toLocalFile()
		self.load_file(fname)

class EnergyLevelsTrendWindow(EQWidget):
	def __init__(self, id, parent=None):
		super().__init__(id)
		self.setWindowTitle("Energy Levels Trend")
		self.setAcceptDrops(True)
		
		self.dataframe = pd.DataFrame(columns=pyckett.egy_dtypes.keys()).astype(pyckett.egy_dtypes)
		self.dataframe_cache = {}
		self.reduced_energy_input = QQ(QLineEdit, "energylevelstrendwindow_reducedenergy", placeholder="Specify formula for reduced energy (depending on egy, qn1, ..., qn6)")
		self.file_label = QQ(QLabel)
		
		
		self.series_selector = None
		class CustomSeriesSelector(SeriesSelector):
			def changed(self, *args, **kwargs):
				super().changed(*args, **kwargs)
				self.parent.series_selector_changed()
		self.series_selector = CustomSeriesSelector(self)
		self.series_selector.layout.setColumnStretch(8, 0)
		
		layout = QVBoxLayout()
		layout.setAlignment(Qt.AlignmentFlag.AlignCenter)
		
		self.setLayout(layout)
		
		tmp_hlayout = QHBoxLayout()
		layout.addLayout(tmp_hlayout)
		
		tmp = QGroupBox("Options")
		tmp_layout = QGridLayout()
		row = 0
		tmp_layout.addWidget(QQ(QLabel, text="File: "), row, 0)
		tmp_layout.addWidget(self.file_label, row, 1)
		tmp_layout.addWidget(QQ(QPushButton, text="Open", change=lambda x: self.load_file()), row, 2)
		row = 1
		tmp_layout.addWidget(QQ(QLabel, text="Filter: "), row, 0)
		tmp_layout.addWidget(QQ(QLineEdit, "energylevelstrendwindow_egyquery"), row, 1, 1, 2)
		row = 2
		tmp_layout.addWidget(QQ(QLabel, text="Reduced Energy Formula: "), row, 0)
		tmp_layout.addWidget(self.reduced_energy_input, row, 1, 1, 2)
		row = 3
		tmp_layout.addWidget(QQ(QPushButton, text="Transfer to LWP", change=self.transfer_to_main_plot), row, 1, 1, 1)
		
		tmp.setLayout(tmp_layout)
		tmp_hlayout.addWidget(tmp, 1)
		
		
		tmp = QGroupBox("Transition")
		tmp_layout2 = QVBoxLayout()
		tmp.setLayout(tmp_layout2)
		tmp_layout2.addWidget(self.series_selector)
		tmp_hlayout.addWidget(tmp)
		
		tmp_layout = QHBoxLayout()
		
		self.trend_plots = [TrendPlot(get_data=self.get_data, get_init_level=lambda i=i: self.get_init_level(i)) for i in ("upper", "lower")]
		for trend_plot in self.trend_plots:
			tmp_layout.addWidget(trend_plot)
		
		layout.addLayout(tmp_layout, 1)
		
		self.series_selector_changed()

	def load_file(self, fname=None):
		if fname is None:
			fname = QFileDialog.getOpenFileName(None, 'Choose Egy File to load',"")[0]
		if fname:
			self.dataframe = pyckett.egy_to_df(fname)
			self.dataframe_cache = {}
			self.file_label.setText(os.path.split(fname)[1])
			
			for trend_plot in self.trend_plots:
				trend_plot.run_fit()

	def get_init_level(self, level_type):
		state = self.series_selector.get_values()
		
		if level_type == "upper":
			return(state["qnus"])
		else:
			return(state["qnls"])

	def get_data(self):
		query = main.config["energylevelstrendwindow_egyquery"].strip()
		
		if self.dataframe_cache.get(query) != query:
			if query:
				filtered_data = self.dataframe.query(query)
			else:
				filtered_data = self.dataframe
			self.dataframe_cache = {"query": query, "data": filtered_data}
		else:
			filtered_data = self.dataframe_cache["data"]
			
		return(filtered_data)
	
	def series_selector_changed(self):
		if self.series_selector is None:
			return
		state = self.series_selector.get_values()
		
		incr_values = state["qndiffs"] if state["increase_freely"] else state["qnincrs"]
		qnus = state["qnus"]
		qnls = state["qnls"]
		
		is_one_series = ( np.array(qnus) == np.array(qnls) + np.array(incr_values) ).all()
		self.trend_plots[1].setHidden(is_one_series)
		
		
		for i, incr in enumerate(incr_values):
			if incr:
				unique_qn = f"qn{i+1}"
				break
		else:
			main.notification("No quantum number is changing. Fit will not be sensible.")
		
		for trend_plot in self.trend_plots:
			trend_plot.unique_qn = unique_qn
			trend_plot.last_level = None
			
			trend_plot.axs[1].set_xlabel(unique_qn)
			trend_plot.run_fit()
		
		if any([x.fit_data[unique_qn].duplicated().any() for x in self.trend_plots]):
			main.notification(f"There are duplicate values in your data for your unique quantum number {unique_qn}. You might want to reset the used levels.")
			return
	
	def transfer_to_main_plot(self):
		state = self.series_selector.get_values()
		
		noq = main.config["series_qns"]
		incr_values = state["qndiffs"] if state["increase_freely"] else state["qnincrs"]
		qnus = state["qnus"]
		qnls = state["qnls"]
		
		incr_values, qnus, qnls = np.array(incr_values[:noq]), np.array(qnus[:noq]), np.array(qnls[:noq])
		
		is_one_series = ( np.array(qnus) == np.array(qnls) + np.array(incr_values) ).all()
		
		i0, i1 = (0, 0) if is_one_series else (0, 1)
		df_upper_series, df_lower_series = self.trend_plots[i0].fit_data.copy(), self.trend_plots[i1].fit_data.copy()
		
		unique_qn = self.trend_plots[i0].unique_qn
		unique_qn_index = int(unique_qn[2:]) - 1
		
		if df_upper_series[unique_qn].duplicated().any() or df_lower_series[unique_qn].duplicated().any():
			main.notification(f"There are duplicate values in your data for your unique quantum number {unique_qn}.")
			return
		
		unique_qn_diff = int(incr_values[unique_qn_index])
		unique_qn_value_reference = qnus[unique_qn_index]
		df_upper_series["index"] = df_upper_series.eval(f"({unique_qn} - {unique_qn_value_reference}) / {unique_qn_diff}").astype(int)
		unique_qn_value_reference = qnls[unique_qn_index]
		df_lower_series["index"] = df_lower_series.eval(f"({unique_qn} - {unique_qn_value_reference}) / {unique_qn_diff}").astype(int)
		
		merge_df = pd.merge(df_upper_series, df_lower_series, how="inner", on="index").sort_values("index")
		merge_df["transition"] = (merge_df["egy_x"] - merge_df["egy_y"]) * main.config["energylevelstrendwindow_egycatconversionfactor"]
		
		positions = merge_df["transition"].tolist()

		if main.config["energylevelstrendwindow_egyqns"]:
			qnu_labels = [f"qn{i+1}_x" for i in range(noq)]
			qnl_labels = [f"qn{i+1}_y" for i in range(noq)]
			qnus = merge_df[qnu_labels].values.tolist()
			qnls = merge_df[qnl_labels].values.tolist()
			quantum_numbers = list(zip(qnus, qnls))
		else:
			indices = merge_df["index"].values
			quantum_numbers = [[(qnus + incr_values * i).tolist(), (qnls + incr_values * i).tolist()] for i in indices]
		
		reference = main.config["series_references"][main.config["series_currenttab"]]
		reference["method"] = "List"
		reference["list"] = {"qns": quantum_numbers, "xs": positions, "i0": 0}

		main.mainwindow.referenceserieswindow.set_state(main.config["series_references"])
		main.plotwidget.set_data()
		main.plotwidget.activateWindow()
	
	def dragEnterEvent(self, event):
		mimeData = event.mimeData()
		if mimeData.hasUrls() and len(mimeData.urls()) == 1:
			event.accept()
		else:
			event.ignore()

	def dropEvent(self, event):
		url = event.mimeData().urls()[0]
		fname = url.toLocalFile()
		self.load_file(fname)

class SpectraResolverWindow(EQWidget):
	def __init__(self, id, parent=None):
		super().__init__(id, parent)
		self.setWindowTitle("Spectra Resolver")
		self.list_ = QListWidget()
		self.list_.setDragDropMode(QAbstractItemView.DragDropMode.InternalMove)
		self.list_.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)

		self.label = QQ(QLabel, text="Ready", wordwrap=True)

		self.apply_order_button = QQ(QPushButton, text="Save overlap free")
		self.apply_order_button.clicked.connect(self.solveoverlap)

		self.layout = QVBoxLayout()
		self.layout.addWidget(self.list_)
		self.layout.addWidget(self.label)
		self.layout.addWidget(QQ(QLabel, wordwrap=True, text="Above all loaded Spectra are shown. Sort them by priority by dragging them with the mouse. Then press 'Save overlap free' to create an overlap free spectrum of all the above spectra."))
		self.layout.addWidget(self.apply_order_button)
		self.setLayout(self.layout)

		self.fill_list()
		main.signalclass.fileschanged.connect(self.fill_list)
		main.signalclass.overlapend.connect(lambda: self.label.setText(f"Ready"))
		main.signalclass.overlapend.connect(lambda: self.apply_order_button.setEnabled(True))
		main.signalclass.overlapindicator.connect(self.label.setText)

	def fill_list(self):
		self.apply_order_button.setEnabled(False)
		main.app.processEvents()
		try:
			entries_pre = [self.list_.item(i).text() for i in range(self.list_.count())]
			entries_now = list(main.config["files_exp"].keys())

			entries_pre = [x for x in entries_pre if x in entries_now]
			entries_now = [x for x in entries_now if x not in entries_pre]

			self.list_.clear()

			for fname in entries_pre + entries_now:
				tmp = QListWidgetItem(fname)
				self.list_.addItem(tmp)
		finally:
			self.apply_order_button.setEnabled(True)

	def solveoverlap(self):
		fname_savefile = QFileDialog.getSaveFileName(None, 'Choose File to Save to',"")[0]
		if fname_savefile == "" or len(main.exp_df) == 0 or len(main.config["files_exp"]) == 0:
			return
		self.apply_order_button.setEnabled(False)
		self.label.setText("Working")
		self.solveoverlap_core(fname_savefile)

	@threading_d
	def solveoverlap_core(self, fname_savefile):
		try:
			df = main.get_visible_data("exp", force_all=True)
			df["keep"] = 1
			i_keep = df.columns.get_loc("keep")

			ranked_files = [self.list_.item(i).text() for i in range(self.list_.count())]
			results = []

			for i, fname in enumerate(ranked_files):
				main.signalclass.overlapindicator.emit(f"Working on file {i+1} of {len(ranked_files)}")
				min, max = main.config["files_exp"][fname]["xrange"]
				i_start, i_stop = df["x"].searchsorted(min, side="left"), df["x"].searchsorted(max, side="right")

				tmp_df = df.iloc[i_start:i_stop]
				df.iloc[i_start:i_stop, i_keep] = 0
				tmp_df = tmp_df[tmp_df["filename"] == fname]
				results.append(tmp_df.copy())
				df = df[df["keep"] == 1]

			main.signalclass.overlapindicator.emit(f"Working on saving the results")
			df = pd.concat(results, ignore_index=True)
			df.sort_values("x", inplace=True, kind="merge")
			df[["x", "y"]].to_csv(fname_savefile, header=None, index=None, sep=chr(main.config["flag_separator"]))
		finally:
			main.signalclass.overlapend.emit()

class ReportWindow(EQWidget):
	def __init__(self, id, parent=None):
		super().__init__(id, parent)
		self.setWindowTitle("Report Window")

		layout = QVBoxLayout()
		self.setLayout(layout)

		layout.addWidget(QQ(QPlainTextEdit, "reportwindow_query", maxHeight=40, placeholder="Query text to filter assignments. Use qnu1, ..., qnu6 and qnl1, ..., qnl6 for the quantum numbers. Other possible values are x, error, weight, comment, and filename."))
		layout.addWidget(QQ(QCheckBox, "reportwindow_blends", text="Blends"))
		layout.addWidget(QQ(QPushButton, text="Create Report", change=self.create_report))
		self.reportfield = QQ(QTextEdit, readonly=True)
		self.reportfield.setFontFamily("Courier")
		layout.addWidget(self.reportfield)

	def create_report(self):
		results = {}
		report = []

		lin_df = main.get_visible_data("lin")
		cat_df = main.get_visible_data("cat")
		noq = main.config["series_qns"]
		qns_visible = [f"qn{ul}{n+1}" for ul in ("u", "l") for n in range(noq)]
		blends = main.config["reportwindow_blends"]

		query = main.config["reportwindow_query"].strip()
		if query:
			lin_df.query(query, inplace=True)

		results["not"] = len(lin_df)
		results["nol"] = lin_df["x"].nunique()
		results["non"] = len(lin_df[lin_df["weight"] == 0])
		results["nod"] = sum(lin_df.duplicated(subset=qns_visible))

		report.append((f"Transitions:", results["not"]))
		report.append((f"Lines:", results["nol"]))
		report.append((f"Unweighted Transitions:", results["non"]))
		report.append((f"Duplicates:", results["nod"]))
		report.append(("", ""))

		results[f"x_min"] = lin_df["x"].min()
		results[f"x_max"] = lin_df["x"].max()
		report.append((f"Frequency min:", results['x_min']))
		report.append((f"Frequency max:", results['x_max']))

		for i in range(noq):
			for ul in ("u", "l"):
				tag = f"qn{ul}{i+1}"
				results[f"{tag}_min"] = lin_df[tag].min()
				results[f"{tag}_max"] = lin_df[tag].max()
				report.append((f"{tag} min:", results[f'{tag}_min']))
				report.append((f"{tag} max:", results[f'{tag}_max']))


		results["nocd"]= sum(cat_df.duplicated(subset=qns_visible))
		if results["nocd"]:
			cat_df = cat_df.drop_duplicates(qns_visible, keep="first")

		df = pd.merge(lin_df, cat_df, how="inner", on=qns_visible)
		df.rename(columns={"x_x": "x_lin", "x_y": "x_cat", "error_x": "error_lin", "error_y": "error_cat", "filename_x": "filename_lin", "filename_y": "filename_cat"}, inplace=True)
		df.reset_index(drop=True, inplace=True)

		if blends:
			mask = df["weight"] != 0
			df_view = df[mask]
			tmp_dict = df_view.groupby(df_view.x_lin).apply(lambda x: np.average(x.x_cat, weights=x.weight)).to_dict()
			df["x_cat"] = df["x_lin"].map(lambda x: tmp_dict.get(x, x))

		df["diff"] = abs(df["x_cat"] - df["x_lin"])
		df["wdiff"] = abs(df["x_cat"] - df["x_lin"])/abs(df["error_lin"])

		results["nom"] = len(df)
		results["max_deviation"] = df["diff"].max()
		results["max_wdeviation"] = df["wdiff"].max()

		report.append(("", ""))
		report.append(("Assignments with matching predict.:", results['nom']))
		report.append(("Assignments w/o  matching predict.:", results['not']-results['nom']))
		report.append(("Highest absolute deviation:", results['max_deviation']))
		report.append(("Highest relative deviation:", results['max_wdeviation']))

		results["rms"] = np.sqrt( (df["diff"]**2).mean() )
		results["wrms"] = np.sqrt( (df["wdiff"]**2).mean() )

		report.append(("", ""))
		report.append(("RMS:", results['rms']))
		report.append(("WRMS:", results['wrms']))

		report = [f"{title: <36}{value:16.4f}" if title else "" for title, value in report]

		if results["nocd"]:
			report.append(f"\nWARNING: {results['nocd']} duplicates were found in your predictions. Each first occurence was kept.")
		if results["not"] != results["nom"]:
			report.append(f"\nWARNING: {results['not']-results['nom']} assignments have no matching prediction. This affects i.a. the RMS and WRMS.")
		if any(df["error_lin"] == 0):
			report.append(f"\nWARNING: Some errors (uncertainties) of your assignments are zero. This leads to infinity values for the relative deviation and the WRMS. Consider using 'error != 0' in the query field.")

		report = "\n".join(report)

		tmp = (self.reportfield.verticalScrollBar().value(), self.reportfield.horizontalScrollBar().value())
		self.reportfield.setText(report)
		self.reportfield.verticalScrollBar().setValue(tmp[0])
		self.reportfield.horizontalScrollBar().setValue(tmp[1])

class FigureWindow(EQWidget):
	def __init__(self, id, parent=None):
		super().__init__(id, parent)
		self.setWindowTitle("Figure Window")

		layout = QHBoxLayout()
		self.setLayout(layout)

		row_index = 0

		tmplayout = QGridLayout()
		layout.addLayout(tmplayout)

		tmplayout.addWidget(QQ(QLabel, text="Size:"), 0, 0)
		tmplayout2 = QHBoxLayout()
		tmplayout2.addWidget(QQ(QDoubleSpinBox, "figurewindow_width", range=(0, None), minWidth=60))
		tmplayout2.addWidget(QQ(QLabel, text=" x "))
		tmplayout2.addWidget(QQ(QDoubleSpinBox, "figurewindow_height", range=(0, None), minWidth=60))
		tmplayout2.addWidget(QQ(QComboBox, "figurewindow_unit", options=("cm", "inch")))
		tmplayout.addLayout(tmplayout2, 0, 1)

		row_index += 1
		tmplayout.addWidget(QQ(QLabel, text="DPI: "), row_index, 0)
		tmplayout.addWidget(QQ(QDoubleSpinBox, "figurewindow_dpi", range=(0, None)), row_index, 1)

		row_index += 1
		tmplayout.addWidget(QQ(QLabel, text="Font-Size: "), row_index, 0)
		tmplayout.addWidget(QQ(QDoubleSpinBox, "figurewindow_fontsize", range=(0, None)), row_index, 1)

		row_index += 1
		tmplayout.addWidget(QQ(QLabel, text="x-Label: "), row_index, 0)
		tmplayout.addWidget(QQ(QLineEdit, "figurewindow_xlabel"), row_index, 1)

		row_index += 1
		tmplayout.addWidget(QQ(QLabel, text="x-Label Padding: "), row_index, 0)
		tmplayout.addWidget(QQ(QDoubleSpinBox, "figurewindow_xlabelpadding", range=(None, None)), row_index, 1)

		row_index += 1
		tmplayout.addWidget(QQ(QLabel, text="y-Label: "), row_index, 0)
		tmplayout.addWidget(QQ(QLineEdit, "figurewindow_ylabel"), row_index, 1)

		row_index += 1
		tmplayout.addWidget(QQ(QLabel, text="y-Label Padding: "), row_index, 0)
		tmplayout.addWidget(QQ(QDoubleSpinBox, "figurewindow_ylabelpadding", range=(None, None)), row_index, 1)

		row_index += 1
		tmplayout.addWidget(QQ(QLabel, text="Annotation: "), row_index, 0)
		tmplayout.addWidget(QQ(QComboBox, "figurewindow_annotation", options=("Like main plot", "Custom")), row_index, 1)

		visible = main.config["figurewindow_annotation"] == "Custom"
		ca_widgets = [QQ(QLabel, text="Custom annotation: ", visible=visible), QQ(QLineEdit, "figurewindow_customannotation", visible=visible)]
		row_index += 1
		tmplayout.addWidget(ca_widgets[0], row_index, 0)
		tmplayout.addWidget(ca_widgets[1], row_index, 1)

		main.config.register("figurewindow_annotation", lambda: ([elem.setVisible(main.config["figurewindow_annotation"] == "Custom") for elem in ca_widgets]))

		row_index += 1
		tmplayout.addWidget(QQ(QLabel, text="Layout: "), row_index, 0)
		tmplayout.addWidget(QQ(QComboBox, "figurewindow_border", options=("Automatic", "Individual")), row_index, 1)

		border_widgets = []
		visible = main.config["figurewindow_border"] == "Individual"
		for side in ("left", "right", "top", "bottom"):
			row_index += 1
			tmp_label = QQ(QLabel, text=f"Border {side}: ", visible=visible)
			tmplayout.addWidget(tmp_label, row_index, 0)
			tmp_input = QQ(QDoubleSpinBox, f"figurewindow_{side}border", range=(0, 1), visible=visible)
			tmplayout.addWidget(tmp_input, row_index, 1)
			border_widgets.extend((tmp_label, tmp_input))

		main.config.register("figurewindow_border", lambda: ([elem.setVisible(main.config["figurewindow_border"] == "Individual") for elem in border_widgets]))

		row_index += 1
		tmplayout.addWidget(QQ(QPushButton, text="Create Figure", change=self.create_figure), row_index, 0, 1, 2)
		row_index += 1
		tmplayout.addWidget(QQ(QPushButton, text="Save Figure", change=self.save_figure), row_index, 0, 1, 2)

		row_index += 1
		tmplayout.setRowStretch(row_index, 1)

		self.fig = None
		self.scene = QGraphicsScene(self)
		self.view = QGraphicsView(self.scene)
		self.view.setStyleSheet("background:transparent;")

		layout.addWidget(self.view, 2)

		self.layout = layout

		QShortcut(
			QKeySequence(Qt.Key.Key_ZoomIn),
			self.view,
			activated=self.zoom_in,
		)

		QShortcut(
			QKeySequence(Qt.Key.Key_ZoomOut),
			self.view,
			activated=self.zoom_out,
		)

	def wheelEvent(self,event):
		steps = event.angleDelta().y() // 120
		tmp_funct = self.zoom_in if steps > 0 else self.zoom_out
		for i in range(abs(steps)):
			tmp_funct()

	def zoom_in(self):
		transform = QTransform()
		transform.scale(1.2, 1.2)

		transform = self.view.transform() * transform
		self.view.setTransform(transform)

	def zoom_out(self):
		transform = QTransform()
		transform.scale(1.2, 1.2)

		transform, invertible = transform.inverted()

		if invertible:
			transform = self.view.transform() * transform
			self.view.setTransform(transform)

	def save_figure(self):
		if not self.fig:
			return
		fname = QFileDialog.getSaveFileName(None, 'Choose file to save figure to')[0]
		if not fname:
			return

		self.fig.savefig(fname)

	def create_figure(self):
		initial_fontsize = matplotlib.rcParams["font.size"]
		try:
			width, height, unit = main.config["figurewindow_width"], main.config["figurewindow_height"], main.config["figurewindow_unit"]
			if unit == "cm":
				width  /= 2.54
				height /= 2.54

			matplotlib.rcParams["font.size"] = main.config["figurewindow_fontsize"]

			rows = max(main.config["plot_rows"], 1)
			cols = max(main.config["plot_cols"], 1)

			self.fig = figure.Figure(figsize=(width, height), dpi=main.config["figurewindow_dpi"], **main.config["figurewindow_figkwargs"])
			self.axs = self.fig.subplots(rows, cols, gridspec_kw=main.config["plot_matplotlibkwargs"], squeeze=False)

			cat_df = main.get_visible_data("cat", scale=False)
			xpos, qns = main.plotwidget.get_positions(return_qns=True, cat_df=cat_df)
			widths = main.plotwidget.get_widths()
			widths += (widths == 0)*10

			offsets = main.plotwidget.get_offset(all=True)
			xmax = xpos + offsets + widths/2
			xmin = xpos + offsets - widths/2

			# set data and set y-ranges of data
			scaling = main.config["plot_yscale"]

			datatypes = ("exp", "cat", "lin")
			dataframes = [main.get_visible_data("exp", scale=False), cat_df, main.get_visible_data("lin", scale=False)]
			files_dicts = [main.config[f"files_{type}"] for type in datatypes]

			minindices = {type: dataframe["x"].searchsorted(xmin, side="left") for type, dataframe in zip(datatypes, dataframes)}
			maxindices = {type: dataframe["x"].searchsorted(xmax, side="right") for type, dataframe in zip(datatypes, dataframes)}

			scalingfactordicts = {type: {file: main.config[f"files_{type}"][file].get("scale", 1) for file in main.config[f"files_{type}"].keys()} for type in ("exp", "cat")}

			ann_dict = main.config["plot_annotations_dict"]


			for i in range(rows):
				for j in range(cols):
					ax = self.axs[i, j]
					ax.yaxis.set_visible(False)
					ax.xaxis.set_visible(False)
					ax.set_xlim(xmin[i, j], xmax[i, j])

					tmp_qns = qns[i, j]
					color = matplotlib.rcParams['text.color']

					if main.config["figurewindow_annotation"] == "Custom":
						reference = main.plotwidget.get_series_reference(j)
						if reference["method"] == "Transition":
							qn_dict = {f"qn{ul}{i+1}": tmp_qns[0 + (ul=="l")][i] for ul in "ul" for i in range(len(tmp_qns[0]))}
						else:
							qn_dict = {}
						text = main.config["figurewindow_customannotation"].format(**{"x": xpos[i, j], "i": i, "j": j, "rows": rows, "cols": cols, **qn_dict})
					elif main.config["series_annotate_xs"] or tmp_qns is None:
						text = f"{{:{main.config['series_annotate_fmt']}}}".format(xpos[i, j])
					else:
						text = f"{', '.join([str(qn) if qn != pyckett.SENTINEL else '-' for qn in tmp_qns[0]])} ← {', '.join([str(qn) if qn != pyckett.SENTINEL else '-' for qn in tmp_qns[1]])}"

						lin_df = main.get_visible_data("lin")
						if len(lin_df.query(" and ".join([f"qnu{i+1} == {qn}" for i, qn in enumerate(tmp_qns[0])] + [f"qnl{i+1} == {qn}" for i, qn in enumerate(tmp_qns[1])]))):
							color = main.config["color_lin"]

					ax.text(**ann_dict, s=text, transform=ax.transAxes, color=color)


					for datatype, dataframe, files in zip(("exp", "cat", "lin"), dataframes, files_dicts):
						minindex, maxindex = minindices[datatype][i][j], maxindices[datatype][i][j]
						dataframe = dataframe.iloc[minindex:maxindex].copy()

						xs = dataframe["x"].to_numpy()
						if datatype == "lin":
							ys = 0*xs
						else:
							ys = dataframe["y"]*dataframe["filename"].replace(scalingfactordicts[datatype]).to_numpy()

						if datatype == "exp":
							if scaling == "Per Plot":
								if len(dataframe):
									yrange_exp = [ys.min(), ys.max()]
								else:
									yrange_exp = [-1, 1]

							filenames = dataframe["filename"].to_numpy()
							unique_filenames = np.unique(filenames)

							for unique_filename in unique_filenames:
								mask = (filenames == unique_filename)
								tmp_xs, tmp_ys = xs[mask], ys[mask]

								ax.plot(tmp_xs, tmp_ys, **{"color": files[unique_filename]["color"], **main.config["figurewindow_expkwargs"]})

						elif datatype == "cat":
							if scaling == "Per Plot":
								if len(dataframe):
									yrange_cat = [dataframe["y"].min(), dataframe["y"].max()]
								else:
									yrange_cat = [-1, 1]
								ys = ys*yrange_exp[1]/yrange_cat[1]
							elif scaling in ["Global", "Custom"]:
								ys = ys*main.config["plot_expcat_factor"]*10**main.config["plot_expcat_exponent"]

							segs = np.array(((xs, xs), (ys*0, ys))).T
							colors = create_colors(dataframe, files, xpos[i, j])
							coll = matplotlib.collections.LineCollection(segs, **{"colors": colors, **main.config["figurewindow_catkwargs"]})
							ax.add_collection(coll)
						elif datatype == "lin":
							colors = create_colors(dataframe, files, lin=True)
							ax.scatter(xs, ys, **{"color": colors, **main.config["figurewindow_linkwargs"]})

					if scaling == "Per Plot":
						yrange = yrange_exp
					elif scaling == "Global":
						yrange = main.yrange_exp
					else:
						yrange = (main.config["plot_yscale_min"], main.config["plot_yscale_max"])

					margin = main.config["plot_ymargin"]

					yrange = [yrange[0]-margin*(yrange[1]-yrange[0]), yrange[1]+margin*(yrange[1]-yrange[0])]
					if np.isnan(yrange[0]) or np.isnan(yrange[1]) or yrange[0] == yrange[1]:
						yrange = [-1,+1]
					ax.set_ylim(yrange)

			# set ticks for bottom row
			i = -1
			for j in range(cols):
				ax = self.axs[i, j]
				ax.xaxis.set_visible(True)
				ticks = np.linspace(xmin[i, j], xmax[i, j], main.config["plot_ticks"])

				reference = main.plotwidget.get_series_reference(j)
				if reference["method"] == "Fortrat":
					ticklabels = symmetric_ticklabels(ticks/2/reference["fortrat"]["jstart"])
				elif main.config["plot_offsetticks"] == 1:
					ticklabels = symmetric_ticklabels(ticks-xpos[i, j])
				elif main.config["plot_offsetticks"] == 0:
					ticklabels = symmetric_ticklabels(ticks)
				else:
					ticklabels = [f"{x:.2e}".replace("e+00", "").rstrip("0").rstrip(".") for x in ticks]

				if j!=0 and len(ticks)>1:
					ticks = ticks[1:]
					ticklabels = ticklabels[1:]
				ax.set_xticks(ticks)
				ax.set_xticklabels(ticklabels)

			xlabel, xlabelpad = main.config["figurewindow_xlabel"], main.config["figurewindow_xlabelpadding"]
			ylabel, ylabelpad = main.config["figurewindow_ylabel"], main.config["figurewindow_ylabelpadding"]
			self.label_ax = shared_labels(self.fig, xlabel, ylabel, xlabelpad, ylabelpad)

			if main.config["figurewindow_border"] == "Automatic":
				self.fig.set_tight_layout(True)
			else:
				self.fig.subplots_adjust(**{side: main.config[f"figurewindow_{side}border"] for side in ("left", "right", "top", "bottom")})

			self.update_figure()
		except Exception as E:
			main.notification(f"<span style='color:#ff0000;'>ERROR</span>: : The figure could not be created. The error message is {E}")
			raise
		finally:
			matplotlib.rcParams["font.size"] = initial_fontsize

	def update_figure(self):
		for item in self.scene.items():
			self.scene.removeItem(item)
		self.plot_canvas = FigureCanvas(self.fig)
		self.scene.addWidget(self.plot_canvas)
		rect = self.scene.itemsBoundingRect()
		self.scene.setSceneRect(rect)

class ConfigWindow(EQWidget):
	def __init__(self, id, parent=None):
		super().__init__(id, parent)
		self.setWindowTitle("Config")

		vbox = QVBoxLayout()
		scrollarea = QScrollArea()
		widget = QWidget()
		layout = QGridLayout()

		self.updating = True

		scrollarea.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
		scrollarea.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
		scrollarea.setWidgetResizable(True)

		tmp_layout = QHBoxLayout()
		tmp_layout.addWidget(QQ(QPushButton, text="Save as default", change=lambda: main.saveoptions()))
		completer = QCompleter(main.config.keys())
		completer.setCaseSensitivity(Qt.CaseSensitivity.CaseInsensitive)
		tmp_layout.addWidget(QQ(QLineEdit, placeholder="Search", completer=completer, change=lambda x: self.search(x)))
		tmp_layout.addStretch(1)

		vbox.addLayout(tmp_layout)
		self.widgets = {}

		i = 1
		for key, value in main.config.items():
			text = json.dumps(value) if isinstance(value, (dict, list, tuple)) else str(value)
			tmp_input = QQ(QLineEdit, value=text, change=lambda text, key=key: self.set_value(key, text))
			tmp_oklab = QQ(QLabel, text="Good")
			tmp_label = QQ(QLabel, text=key)

			self.widgets[key] = (tmp_input, tmp_oklab, tmp_label)
			layout.addWidget(tmp_label, i+1, 0)
			layout.addWidget(tmp_input, i+1, 1)
			layout.addWidget(tmp_oklab, i+1, 2)
			i += 1

		layout.setRowStretch(i+1, 1)

		widget.setLayout(layout)
		scrollarea.setWidget(widget)
		vbox.addWidget(scrollarea)

		self.setLayout(vbox)

		self.updating = False

	def show(self, *args, **kwargs):
		self.timer = QTimer(self)
		self.timer.timeout.connect(self.get_values)
		self.timer.start(200)
		return(super().show(*args, **kwargs))

	def search(self, text):
		for key, value in self.widgets.items():
			if text.lower() in key or text.lower() in value[0].text():
				hidden = False
			else:
				hidden = True
			value[0].setHidden(hidden)
			value[1].setHidden(hidden)
			value[2].setHidden(hidden)

	def get_values(self):
		self.updating = True
		for key, (input, oklabel, label) in self.widgets.items():
			value = main.config[key]
			if input.hasFocus() or self.widgets[key][1].text() == "Bad":
				continue
			if isinstance(value, (dict, list, tuple)):
				input.setText(json.dumps(value))
			else:
				input.setText(str(value))
		self.updating = False

	def set_value(self, key, value):
		if self.updating:
			return
		converter = config_specs.get(key)
		if converter:
			converter = converter[1]
		input, oklab, label = self.widgets[key]

		try:
			if converter is None:
				pass
			elif converter in (dict, list, tuple):
				value = json.loads(value)
			elif converter == bool:
				value = True if value in ["True", "1"] else False
			else:
				value = converter(value)
			main.config[key] = value
			oklab.setText("Good")
		except Exception as E:
			oklab.setText("Bad")


	def closeEvent(self, *args, **kwargs):
		self.timer.stop()
		return super().closeEvent(*args, **kwargs)

class CalibrateSpectrumWindow(EQWidget):
	def __init__(self, id, parent=None):
		super().__init__(id, parent)
		self.setWindowTitle("Spectrum Calibration")
		
		vbox = QVBoxLayout()
		self.setLayout(vbox)
		
		self.cal_df = None
		self.calibration_points = []

		self.fig = figure.Figure(dpi=main.config["plot_dpi"])
		self.plotcanvas = FigureCanvas(self.fig)
		self.plotcanvas.setMinimumHeight(200)
		self.plotcanvas.setMinimumWidth(200)
		
		self.axs = self.fig.subplots(3, sharex=True, gridspec_kw = {'wspace':0, 'hspace':0})
		self.axs[0].get_xaxis().set_visible(False)
		self.axs[1].get_xaxis().set_visible(False)
		self.factors  = self.axs[0].scatter([], [], color=main.config["calibratewindow_color"])
		self.exp_line = self.axs[1].plot([], [], color=main.config["color_exp"])[0]
		self.cal_line = self.axs[2].plot([], [], color=main.config["calibratewindow_color"])[0]
	
		self.mpltoolbar = NavigationToolbar2QT(self.plotcanvas, self)

		vbox.addWidget(self.plotcanvas)
		vbox.addWidget(self.mpltoolbar)
		
		vbox.addWidget(QQ(QPlainTextEdit, "calibratewindow_queryexp", maxHeight=60, placeholder="Query text to filter experimental dataframe."))
		vbox.addWidget(QQ(QPlainTextEdit, "calibratewindow_querycat", maxHeight=60, placeholder="Query text to filter predictions dataframe."))

		tmplayout = QGridLayout()
		row_id = 0
		tmplayout.addWidget(QQ(QLabel, text="Function: "), row_id, 0)
		tmplayout.addWidget(QQ(QComboBox, "calibratewindow_lineshape", items=("Gauss", "Lorentz", "Voigt"), minWidth=120), row_id, 1)

		tmplayout.addWidget(QQ(QLabel, text="Color: "), row_id, 2)
		tmplayout.addWidget(QQ(QLineEdit, "calibratewindow_color", minWidth=120, ), row_id, 3)
		
		row_id += 1

		tmplayout.addWidget(QQ(QLabel, text="Derivative: "), row_id, 0)
		tmplayout.addWidget(QQ(QSpinBox, "calibratewindow_derivative", range=(0, 2), minWidth=120), row_id, 1)
		
		row_id += 1

		tmplayout.addWidget(QQ(QLabel, text="Blendwidth: "), row_id, 0)
		tmplayout.addWidget(QQ(QDoubleSpinBox, "calibratewindow_blendwidth", range=(0, None), minWidth=120), row_id, 1)
		
		tmplayout.addWidget(QQ(QLabel, text="Fitwidth: "), row_id, 2)
		tmplayout.addWidget(QQ(QDoubleSpinBox, "calibratewindow_fitwidth", range=(0, None), minWidth=120), row_id, 3)
		
		tmplayout.setColumnStretch(5, 2)
		
		vbox.addLayout(tmplayout)
		
		self.run_button = QQ(QPushButton, text="Run", change=lambda x: self.run())
		self.save_button = QQ(QPushButton, text="Save", change=lambda x: self.save())
		vbox.addWidget(self.run_button)
		vbox.addWidget(self.save_button)
		
		main.signalclass.calibrationstart.connect(lambda: self.run_button.setEnabled(False))
		main.signalclass.calibrationend.connect(lambda: self.after_calibration())
	
	def update_plot(self):
		if len(self.calibration_points):
			self.factors.set_offsets(self.calibration_points)
			self.factors.set_color(main.config["calibratewindow_color"])
			self.axs[0].set_ylim(self.calibration_points[:, 1].min(), self.calibration_points[:, 1].max())
		
		if isinstance(self.cal_df, pd.DataFrame):
			xs = self.cal_df["x"].to_numpy()
			ys_old = self.cal_df["y_old"].to_numpy()
			ys = self.cal_df["y"].to_numpy()
			self.exp_line.set_data(xs, ys_old)
			self.axs[1].set_ylim(ys_old.min(), ys_old.max())
			self.cal_line.set_data(xs, ys)
			self.axs[2].set_ylim(ys.min(), ys.max())
		
		self.axs[0].set_xlim(xs.min(), xs.max())
		self.plotcanvas.draw()
	
	
	@threading_d
	@synchronized_d(locks["calibration"])
	def run(self):
		main.signalclass.calibrationstart.emit()
		
		try:
			exp_df = main.get_visible_data("exp", scale=False)
			cat_df = main.get_visible_data("cat", scale=False)
			
			if main.config["calibratewindow_queryexp"].strip():
				exp_df = exp_df.query(main.config["calibratewindow_queryexp"])
			if main.config["calibratewindow_querycat"].strip():
				cat_df = cat_df.query(main.config["calibratewindow_querycat"])
			
			lineshape_type = main.config["calibratewindow_lineshape"]
			blendwidth = main.config["calibratewindow_blendwidth"] / 2
			fitwidth = main.config["calibratewindow_fitwidth"] / 2
			fitfunction = lambda *args: lineshape(lineshape_type, main.config["calibratewindow_derivative"], *args[:-1]) + args[-1]
			
			xs = cat_df["x"].to_numpy()
			
			i_starts_exp, i_stops_exp = exp_df["x"].searchsorted(xs - fitwidth, side="left"), exp_df["x"].searchsorted(xs + fitwidth, side="right")
			i_starts_cat, i_stops_cat = cat_df["x"].searchsorted(xs - blendwidth, side="left"), cat_df["x"].searchsorted(xs + blendwidth, side="right")
			
			calibration_points = []
			for i, x0 in enumerate(xs):
				# Sum intensity of all blended peaks
				blends_df = cat_df.iloc[i_starts_cat[i]: i_stops_cat[i]]
				y0 = blends_df["y"].sum()
				
				# Get experimental data for fit
				expfit_df = exp_df.iloc[i_starts_exp[i]: i_stops_exp[i]]
				if len(expfit_df) < 10:
					continue
				
				xs, ys = expfit_df["x"].to_numpy(), expfit_df["y"].to_numpy()
				
				# Fit peaks and get calibration factors for these points
				ymax = ys.max()
				p0 = [x0, ymax, fitwidth/10, fitwidth/10, 0]
				bounds = [
					[x0-0.2, 0, 0, 0, 0],
					[x0+0.2, 2*ymax, fitwidth/2.5, fitwidth/2.5, ymax/2],
				]
				
				if lineshape_type != "Voigt":
					del p0[2]
					del bounds[0][2]
					del bounds[1][2]
				
				
				popt, pcov = optimize.curve_fit(fitfunction, xs, ys, p0=p0, bounds=bounds)
				calibration_points.append((popt[0], y0/popt[1]))
				
				# @Luis: Perform checks to see if fit was sensible
					
			if not len(calibration_points):
				calibration_points = np.array([[], []])
				main.notification("There are no calibration points.")
				return
				
			calibration_points = np.array(calibration_points)
			calibration_points = calibration_points[np.argsort(calibration_points[:, 0])]
			
			if not len(exp_df):
				main.notification("There are no experimental points.")
				return
			
			factors = np.interp(exp_df["x"].to_numpy(), calibration_points[:, 0], calibration_points[:, 1], left=None, right=None)
			exp_df["y_old"] = exp_df["y"]
			exp_df["y"] *= factors
			exp_df["y"] /= exp_df["y"].max()
			exp_df["factors"] = factors
			
			main.notification(f"Used {len(calibration_points)} of {len(cat_df)} predictions for calibration.")
			
			self.calibration_points = calibration_points
			self.cal_df = exp_df
		except Exception as E:
			raise
		finally:
			main.signalclass.calibrationend.emit()
	
	def after_calibration(self):
		self.run_button.setEnabled(True)
		self.update_plot()

	def save(self):
		self.save_button.setEnabled(False)
		fname = QFileDialog.getSaveFileName(None, 'Choose file to save spectrum to',"","CSV Files (*.csv);;All Files (*)")[0]
		if fname:
			self.cal_df[["x", "y"]].to_csv(fname, index=None, header=None, sep="\t")
		self.save_button.setEnabled(True)


##
## Dialogs
##
class QNsDialog(QDialog):
	def __init__(self, frequency, df):
		super().__init__()
		QShortcut("Esc", self).activated.connect(lambda: self.predone(0))

		noq = main.config["series_qns"]
		qns = [f"qnu{i+1}" for i in range(6)]+[f"qnl{i+1}" for i in range(6)]
		qns_visible = [f"qn{ul}{n+1}" for ul in ("u", "l") for n in range(noq)]
		self.res = {key: pyckett.SENTINEL for key in qns}

		self.setWindowTitle(f"Choose QNs for transition at {frequency}")
		self.resize(main.config["qnsdialog_width"], main.config["qnsdialog_height"])
		self.frequency = frequency

		layout = QVBoxLayout()
		self.setLayout(layout)
		layout.addWidget(QQ(QLabel, wordwrap=True, text="Enter the quantum numbers in the input fields or choose one of the transitions from the table."))
		layout.addSpacing(10)

		self.sbs = {}
		qnslayout = QGridLayout()
		for i in range(noq):
			widget = QQ(QSpinBox, value=0, range=(0, None))
			qnslayout.addWidget(widget, 1, i)
			self.sbs[f"qnu{i+1}"] = widget

		for i in range(noq):
			widget = QQ(QSpinBox, value=0, range=(0, None))
			qnslayout.addWidget(widget, 2, i)
			self.sbs[f"qnl{i+1}"] = widget

		for i in range(noq):
			lab = QLabel(f"QN {i+1}")
			qnslayout.addWidget(QQ(QLabel, text=f"QN{i+1}"), 0, i)

		qnslayout.setColumnStretch(7, 1)
		layout.addLayout(qnslayout)

		tmp = ["dist", "x", "log y"]
		cols = tmp + qns
		table = QTableWidget()
		table.setColumnCount(len(cols)+1)
		table.setHorizontalHeaderLabels(["Assign"] + cols)
		table.setRowCount(0)
		table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)

		tmp_df = df.copy()
		tmp_df["dist"] = tmp_df["x"] - frequency
		tmp_df["log y"] = np.log10(tmp_df["y"])
		tmp_df["absdist"] = abs(tmp_df["dist"])
		tmp_df.sort_values(by=["absdist"], inplace=True)
		tmp_df.reset_index(drop=True, inplace=True)

		for i, row in tmp_df.head(100).iterrows():
			currRowCount = table.rowCount()
			table.insertRow(currRowCount)
			for j, col in enumerate(cols):
				val = f'{row[col]:{main.config["flag_xformatfloat"]}}'.rstrip("0").rstrip(".")
				table.setItem(currRowCount, j+1, QTableWidgetItem(val))
			tmpd = {key: row[key] for key in qns_visible}
			tmpd["xpre"] = row["x"]
			table.setCellWidget(currRowCount, 0, QQ(QPushButton, text="Assign", change=lambda x, tmpd=tmpd: self.table_save(tmpd)))

		
		for i in range(6):
			table.setColumnHidden(i+ len(tmp) + 1, i>=noq)
			table.setColumnHidden(i+ len(tmp) + 7, i>=noq)

		table.resizeColumnsToContents()
		layout.addWidget(table)

		buttons = QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
		buttonBox = QDialogButtonBox(buttons)
		buttonBox.accepted.connect(lambda: self.selector_save(1))
		buttonBox.rejected.connect(lambda: self.predone(0))
		layout.addWidget(buttonBox)

	def table_save(self, tmpd):
		self.res.update(tmpd)
		self.predone(1)

	def selector_save(self, val):
		dict_ = {key: sb.value() for key, sb in self.sbs.items()}
		dict_["xpre"] = self.frequency
		self.res.update(dict_)
		self.predone(val)

	def save(self):
		return(self.res)

	def predone(self, val):
		main.config["qnsdialog_width"] =	self.geometry().width()
		main.config["qnsdialog_height"] =	self.geometry().height()
		self.done(val)

class ConsoleDialog(QDialog):
	def __init__(self):
		super().__init__()
		QShortcut("Esc", self).activated.connect(lambda: self.predone(0))

		self.setWindowTitle(f"Command Line Dialog")
		self.resize(main.config["commandlinedialog_width"], main.config["commandlinedialog_height"])

		self.tabs = QTabWidget()
		self.tabs.setTabsClosable(True)
		self.tabs.setMovable(True)
		self.tabs.setDocumentMode(True)

		initial_values = main.config["commandlinedialog_commands"]
		if initial_values:
			for title, command in initial_values:
				self.add_tab(title, command)
		else:
			self.add_tab()

		self.tabs.tabCloseRequested.connect(self.close_tab)
		self.tabs.tabBarDoubleClicked.connect(self.renameoradd_tab)
		self.tabs.setCurrentIndex(main.config["commandlinedialog_current"])

		layout = QVBoxLayout()
		self.setLayout(layout)

		layout.addWidget(self.tabs)
		buttons_layout = QHBoxLayout()
		buttons_layout.addStretch()
		buttons_layout.addWidget(QQ(QPushButton, text="Run", change=lambda x: self.predone(1), shortcut="Ctrl+Return"))
		buttons_layout.addWidget(QQ(QPushButton, text="Cancel", change=lambda x: self.predone(0), shortcut="Esc"))
		buttons_layout.addStretch()
		layout.addLayout(buttons_layout)

	def add_tab(self, title="Command", command=""):
		textarea = QQ(QPlainTextEdit, value=command)
		cursor = textarea.textCursor()
		cursor.movePosition(QTextCursor.MoveOperation.End)
		textarea.setTextCursor(cursor)
		self.tabs.addTab(textarea, title)

	def close_tab(self, index):
		tab = self.tabs.widget(index)
		tab.deleteLater()
		self.tabs.removeTab(index)
		if self.tabs.count() == 0:
			self.add_tab()

	def renameoradd_tab(self, index):
		if index == -1:
			self.add_tab()
		elif self.tabs.widget(index) != 0:
			text, ok = QInputDialog().getText(self, "Tab Name","Enter the Tabs Name:")
			if ok and text:
				self.tabs.setTabText(index, text)

	def predone(self, val):
		commands = []
		for i in range(self.tabs.count()):
			tab = self.tabs.widget(i)
			title = self.tabs.tabText(i)
			command = tab.toPlainText()
			commands.append((title, command))

		main.config["commandlinedialog_commands"] = commands
		main.config["commandlinedialog_current"] = self.tabs.currentIndex()

		main.config["commandlinedialog_width"] =	self.geometry().width()
		main.config["commandlinedialog_height"] =	self.geometry().height()
		self.done(val)




##
## Reference and Series Selector
##
class ReferenceSelector(QTabWidget):
	def __init__(self, initial_values, parent=None):
		super().__init__(parent)

		self.state = CDict(lambda: parent.changed(), {
			"method":		"Transition",
			"transition":	None,
			"list": CDict(lambda: parent.changed(), {
			  "qns":		None,
			  "xs":			[],
			  "i0":			0,
			}),
			"expression": CDict(lambda: parent.changed(), {
			  "expression":	"",
			  "N0":			0,
			}),
			"fortrat": CDict(lambda: parent.changed(), {
			  "jstart": 	1,
			  "center": 	8000,
			}),
		})

		if initial_values:
			self.state.update(initial_values)

		plotwidgets = {}
		self.methods = ("Transition", "List", "Expression", "Fortrat")
		for label in self.methods:
			plotwidgets[label] = QWidget()
			self.addTab(plotwidgets[label], label)

		self.setCurrentIndex(self.methods.index(self.state["method"]))
		self.currentChanged.connect(lambda x: self.state.__setitem__("method", self.methods[x]))

		# Tab 1: Transition
		layout = QVBoxLayout()

		text = self.state["transition"].get("file") if self.state.get("transition") else None
		text = os.path.split(text)[-1] if text else "All"
		tooltip = text if text else "The first transition out of all visible files will be used"
		self.activefilebutton = QQ(QPushButton, text=text, tooltip=tooltip, change=lambda x: self.series_selector.select_file())
		tmplayout = QHBoxLayout()
		tmplayout.addWidget(QQ(QLabel, text="Active File: "))
		tmplayout.addWidget(self.activefilebutton, 1)


		class CustomSeriesSelector(SeriesSelector):
			def select_file(self):
				options = ["All"] + list(main.config["files_cat"].keys())
				if self.state.get("file") in options:
					current = options.index(self.state["file"])
				else:
					current = 0
				item, ok = QInputDialog.getItem(self, "Choose File", f"Limit this series to the following file:", options, current=current, editable=False)
				if not (ok and item):
					return

				if item == "All":
					item = None
					itemshort = "All"
				else:
					itemshort = os.path.split(item)[-1]

				self.state["file"] = item
				self.parent.activefilebutton.setText(itemshort)
				self.parent.activefilebutton.setToolTip(item if item else "The first transition out of all visible files will be used")
				self.changed()

			def changed(self):
				super().changed()
				self.parent.state["transition"] = self.get_values()

		self.series_selector = CustomSeriesSelector(self, self.state["transition"])

		layout.addLayout(tmplayout)
		layout.addWidget(self.series_selector)

		button_apply = QQ(QPushButton, text="Apply", change=lambda x: main.plotwidget.reset_offsets())

		checkbox_annotatefrequencies = QQ(QCheckBox, "series_annotate_xs", text="Toggle QNs/Freq")
		label_blend = QLabel("Blend Width: ")
		width_blend = QQ(QDoubleSpinBox, "series_blendwidth", minWidth=85, range=(0, None))
		checkbox_blend = QQ(QCheckBox, "series_blenddialog", text="Blend Dialog", changes=(width_blend.setEnabled, label_blend.setEnabled))
		checkbox_blend.stateChanged.emit(main.config["series_blenddialog"])

		button_box = QHBoxLayout()
		[button_box.addWidget(checkbox_annotatefrequencies), button_box.addWidget(checkbox_blend), button_box.addWidget(label_blend), button_box.addWidget(width_blend), button_box.addStretch(1)]

		layout.addLayout(button_box)
		layout.addWidget(button_apply)
		layout.addStretch(1)
		plotwidgets["Transition"].setLayout(layout)

		# Tab 2: List
		layout = QVBoxLayout()

		self.xsTable = QQ(QTableWidget, rowCount=0, columnCount=2, move=(0, 0))
		self.xsTable.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
		self.xsTable.setHorizontalHeaderLabels(["#", "Frequency"])
		layout.addWidget(self.xsTable)
		if len(self.state["list"]["xs"]):
			self.load_xs_list(values=self.state["list"]["xs"])

		button_open = QQ(QToolButton, text="Open List", change=lambda x: self.load_xs_list(temp=False))
		button_write = QQ(QToolButton, text="Write List", change=lambda x: self.load_xs_list(temp=True))

		label_startat = QLabel("Start at Index: ")
		spinbox_startat = QQ(QSpinBox, value=self.state["list"]["i0"], range=(0, None), singlestep=1, change=lambda: self.state["list"].__setitem__("i0", spinbox_startat.value()))
		self.spinbox_startat = spinbox_startat
		button_apply = QQ(QPushButton, text="Apply", change=lambda x: main.plotwidget.reset_offsets())

		hbox = QHBoxLayout()
		[hbox.addWidget(button_open), hbox.addWidget(button_write), hbox.addStretch(1), hbox.addWidget(label_startat), hbox.addWidget(spinbox_startat), hbox.addWidget(button_apply)]

		layout.addLayout(hbox)

		plotwidgets["List"].setLayout(layout)

		# Tab 3: Expression
		layout = QVBoxLayout()

		self.input_expr = QQ(QPlainTextEdit, value=self.state["expression"]["expression"], placeholder="Enter expression here\ne.g. for a linear molecule with B=4000 use: (N+N0)*4000*2", change=lambda: self.state["expression"].__setitem__("expression", self.input_expr.toPlainText()))
		button_apply = QQ(QPushButton, text="Apply", change=lambda x: main.plotwidget.reset_offsets())

		label_N0 = QQ(QLabel, text="N0: ")
		self.input_N0 = QQ(QSpinBox, value=self.state["expression"]["N0"], range=(None, None), change=lambda: self.state["expression"].__setitem__("N0", self.input_N0.value()))

		layout.addWidget(self.input_expr)
		hbox = QHBoxLayout()
		[hbox.addWidget(label_N0), hbox.addWidget(self.input_N0), hbox.addStretch(1), hbox.addWidget(button_apply)]
		layout.addLayout(hbox)

		plotwidgets["Expression"].setLayout(layout)

		# Tab 4: Fortrat
		layout = QVBoxLayout()

		layout.addWidget(QQ(QLabel, text="J of lowest Plot: "))
		self.jstart = QQ(QSpinBox, value=self.state["fortrat"]["jstart"], singlestep=1, range=(1, None), change=lambda: self.state["fortrat"].__setitem__("jstart", self.jstart.value()))
		layout.addWidget(self.jstart)
		layout.addStretch(1)
		layout.addWidget(QQ(QPushButton, text="Apply", change=lambda x: main.plotwidget.reset_offsets()))

		plotwidgets["Fortrat"].setLayout(layout)

	def get_values(self):
		return(self.state)

	def changed(self):
		pass

	def update(self, newstate):
		self.state.update(newstate)
		self.series_selector.state = self.state["transition"]
		self.series_selector.set_state()

	def load_xs_list(self, values=None, temp=False):
		if not values is None:
			xs = values
			self.state["list"]["xs"] = values

		elif temp:
			line, ok = QInputDialog().getMultiLineText(self, "Specify Custom List",
			"Write list here (delimiters are all whitespace characters, comma and semicolon):")
			if not ok or not line:
				return

			xs = []
			tmp_xs = re.split(r'; |, |\s', line)
			for x in tmp_xs:
				if not x.strip():
					continue
				try:
					xs.append(float(x))
				except ValueError:
					main.notification(f"<span style='color:#eda711;'>WARNING</span>: Could not convert the string '{x}' to a numerical value.")

			self.state["list"].update({
				"qns":		None,
				"xs":		xs,
			})
			self.spinbox_startat.setValue(0)

		else:
			fnames = QFileDialog.getOpenFileNames(self, 'Open Positions List(s)',)[0]
			if len(fnames)==0:
				return
			xs = []
			for fname in fnames:
				with open(fname, "r", encoding="utf-8") as file:
					for line in file:
						line = line.strip()
						if line == "" or line.startswith("#"):
							continue

						tmp = re.split(r'; |, |\s', line)
						for x in tmp:
							try:
								xs.append(float(x))
							except ValueError:
								main.notification(f"<span style='color:#eda711;'>WARNING</span>: Could not convert the string '{x}' to a numerical value.")

			self.state["list"].update({
				"qns":		None,
				"xs":		xs,
			})
			self.spinbox_startat.setValue(0)


		table = self.xsTable
		table.setRowCount(0)
		i=0
		for i, x in enumerate(xs):
			currRowCount = table.rowCount()
			table.insertRow(currRowCount)
			table.setItem(currRowCount, 0, QTableWidgetItem(f"{i}"))
			table.setItem(currRowCount, 1, QTableWidgetItem(f"{x:{main.config['flag_xformatfloat']}}"))
		self.changed()

class SeriesSelector(QWidget):
	def __init__(self, parent, initial_values={}):
		super().__init__(parent)
		self.updating = False
		self.state = {
			"increase_freely":	False,
			"qnus":				[1, 0, 1, 0, 0, 0],
			"qnls":				[0]*6,
			"qnincrs":			[True, False, True, False, False, False],
			"qndiffs":			[1, 0, 1, 0, 0, 0],
		}

		self.parent = parent
		if initial_values:
			self.state.update(initial_values)

		layout = QGridLayout()

		self.labels  = [QQ(QLabel, text=f"QN {x+1}") for x in range(6)]
		self.qnus    = [QQ(QSpinBox, minWidth=40, maxWidth=60, range=(None, None), singlestep=1, change=lambda x: self.changed()) for x in range(6)]
		self.qnls    = [QQ(QSpinBox, minWidth=40, maxWidth=60, range=(None, None), singlestep=1, change=lambda x: self.changed()) for x in range(6)]
		self.qnincrs = [QQ(QCheckBox, text="Incr", minWidth=40, maxWidth=60, change=lambda x: self.changed()) for x in range(6)]
		self.qndiffs = [QQ(QSpinBox, minWidth=40, maxWidth=60, range=(None, None), singlestep=1, change=lambda x: self.changed()) for x in range(6)]

		self.incqns = QQ(QPushButton, text="Increase", change=lambda x: self.incdecqns(+1))
		self.decqns = QQ(QPushButton, text="Decrease", change=lambda x: self.incdecqns(-1))

		self.incnoq = QQ(QPushButton, text="QN+", change=lambda x: self.alternoq(+1))
		self.decnoq = QQ(QPushButton, text="QN-", change=lambda x: self.alternoq(-1))

		self.togglediff = QQ(QToolButton, text="⇆", change=lambda x: self.change_incr_mode())

		for i, widget in enumerate(self.labels + self.qnus + self.qnls):
			layout.addWidget(widget, i//6, i%6)

		for i, cb, diff in zip(range(6), self.qnincrs, self.qndiffs):
			tmp = QHBoxLayout()
			tmp.addWidget(cb)
			tmp.addWidget(diff)
			layout.addLayout(tmp, 4, i)

		for i in range(6):
			layout.setColumnStretch(i, 100)
		layout.addWidget(self.togglediff, 4, 6, 1, 1)

		layout.addWidget(self.incqns, 1, 6, 1, 2)
		layout.addWidget(self.decqns, 2, 6, 1, 2)

		layout.addWidget(self.incnoq, 0, 6)
		layout.addWidget(self.decnoq, 0, 7)

		layout.setRowStretch(6, 10)
		layout.setColumnStretch(8, 1)

		self.layout = layout
		self.setLayout(layout)
		self.set_state()

		main.config.register_widget("series_qns", self.togglediff, lambda: self.alternoq(fromconfig=True))

	def set_state(self):
		self.updating = True
		state = self.state
		visibles = self.get_visible()
		for label, visible in zip(self.labels, visibles):
			label.setVisible(visible)
		for qnu, value, visible in zip(self.qnus, state["qnus"], visibles):
			qnu.setValue(value)
			qnu.setVisible(visible)
		for qnl, value, visible in zip(self.qnls, state["qnls"], visibles):
			qnl.setValue(value)
			qnl.setVisible(visible)
		for qnincr, value, visible in zip(self.qnincrs, state["qnincrs"], visibles):
			qnincr.setChecked(value)
			qnincr.setVisible(visible and not state["increase_freely"])
		for qndiff, value, visible in zip(self.qndiffs, state["qndiffs"], visibles):
			qndiff.setValue(value)
			qndiff.setVisible(visible and state["increase_freely"])
		self.updating = False
		self.changed()

	def change_incr_mode(self):
		self.state["increase_freely"] = not self.state["increase_freely"]
		self.set_state()
		self.changed()

	def incdecqns(self, dir):
		incr_values = (x.value() for x in self.qndiffs) if self.state["increase_freely"] else (x.isChecked() for x in self.qnincrs)

		for qnu, qnl, incr in zip(self.qnus, self.qnls, incr_values):
			qnu.setValue(qnu.value()+dir*incr)
			qnl.setValue(qnl.value()+dir*incr)

	def alternoq(self, value=0, fromconfig=False):
		if value:
			tmp = main.config["series_qns"] + value
		else:
			tmp = main.config["series_qns"]

		if 1 <= tmp <= 6:
			if not fromconfig:
				main.config["series_qns"] = tmp
			self.set_state()
			self.changed()

	def get_visible(self):
		noq = main.config["series_qns"]
		return([True]*noq+[False]*(6-noq))

	def get_values(self):
		return(self.state)

	def changed(self):
		if not self.updating:
			self.state["qnus"] = [x.value() for x in self.qnus]
			self.state["qnls"] = [x.value() for x in self.qnls]
			self.state["qnincrs"] = [x.isChecked() for x in self.qnincrs]
			self.state["qndiffs"] = [x.value() for x in self.qndiffs]

	def wheelEvent(self, event):
		steps = event.angleDelta().y() // 120
		self.incdecqns(steps)




##
## Miscellaneous Classes
##
class FigureCanvas(FigureCanvas):
	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		self.wheelEvent = lambda event: event.ignore()
		self.setStyleSheet('background-color: #00000000')

class ProtPlot(QWidget):
	def __init__(self, parent=None, onrange=False):
		super().__init__(parent)

		self.parent = parent
		self.shortcuts()

		self.i = main.plotwidget.get_current_plot()

		self.center = 5
		self.width = 10

		self.from_current_plot(update=False)
		self.gui()

		self.onrange = onrange
		if onrange:
			self.span_selector = matplotlib.widgets.SpanSelector(self.ax, lambda vmax, vmin: self.on_range(vmax, vmin, self.i), 'horizontal', useblit=True)

	def gui(self):
		layout = QVBoxLayout()
		self.setLayout(layout)

		tmp_layout = QHBoxLayout()
		layout.addLayout(tmp_layout)

		for dir in ("in", "out", "left", "right"):
			tmp_layout.addWidget(QQ(QPushButton, text=dir, change=lambda x, dir=dir: self.move_plot(dir)))

		tmp_layout.addStretch()
		tmp_layout.addWidget(QQ(QPushButton, text="From Current Plot", change=lambda x: self.from_current_plot()))

		self.fig = figure.Figure(dpi=main.config["plot_dpi"])
		self.plot_canvas = FigureCanvas(self.fig)
		layout.addWidget(self.plot_canvas, 6)
		layout.addStretch()

		self.ax = self.fig.subplots()
		self.ax.ticklabel_format(useOffset=False)
		self.exp_line = matplotlib.collections.LineCollection(np.zeros(shape=(0,2,2)), color=main.config["color_exp"])
		self.ax.add_collection(self.exp_line)

	def from_current_plot(self, update=True):
		self.i = main.plotwidget.get_current_plot()
		xrange = main.plotwidget.axs["ax"][self.i].get_xlim()
		self.center = sum(xrange)/2
		if main.config["flag_protplotautowidth"]:
			self.width = xrange[1] - xrange[0]

		if update:
			self.update_plot()

	@synchronized_d(locks["exp_df"])
	def update_plot(self):
		exp_df = main.get_visible_data("exp", xrange=(self.center-self.width/2, self.center+self.width/2), binning=True)

		self.exp_xs = xs = exp_df["x"].to_numpy()
		self.exp_ys = ys = exp_df["y"].to_numpy()

		files = main.config[f"files_exp"]
		if main.config["plot_expasstickspectrum"]:
			segs = np.array(((xs, xs), (ys*0, ys))).T
			colors = create_colors(exp_df, files)
		else:
			# @Luis: optimize this
			filenames = exp_df["filename"].to_numpy()
			unique_filenames = np.unique(filenames)

			segs = []
			colors = []
			for unique_filename in unique_filenames:
				mask = (filenames == unique_filename)
				tmp_xs, tmp_ys = xs[mask], ys[mask]

				segs.append(np.array(((tmp_xs[:-1], tmp_xs[1:]), (tmp_ys[:-1], tmp_ys[1:]))).T)
				colors.extend([files[unique_filename]["color"]]*sum(mask))

			if segs:
				segs = np.concatenate(segs)

		self.exp_line.set(segments=segs, colors=colors)
		self.ax.set_xlim([self.center-self.width/2, self.center+self.width/2])

		yrange = [-1, 1]
		if len(self.exp_ys):
			ymin, ymax = np.min(self.exp_ys), np.max(self.exp_ys)
			width = ymax-ymin
			margin = main.config["plot_ymargin"]*width if width else 1
			yrange = (ymin-margin, ymax+margin)
		self.ax.set_ylim(yrange)

		self.plot_canvas.draw()

	def move_plot(self, dir, factor=None):
		if dir == "in":
			self.width /= 2
		elif dir == "out":
			self.width *= 2
		elif dir == "left":
			self.center -= self.width/2
		elif dir == "right":
			self.center += self.width/2
		elif dir == "sin":
			self.width *= 3/4
		elif dir == "sout":
			self.width /= 3/4
		elif dir == "sleft":
			self.center -= self.width/4
		elif dir == "sright":
			self.center += self.width/4
		elif dir == "wheel" and factor != None:
			self.width *= factor

		self.update_plot()

	def alter_plot(self, pos, dir):
		new_i = list(self.i)
		new_i[pos] += dir
		
		if main.plotwidget.is_valid_plot(new_i):
			self.i = tuple(new_i)
			xrange = main.plotwidget.axs["ax"][self.i].get_xlim()
			self.center = sum(xrange)/2
			self.update_plot()

	def wheelEvent(self,event):
		steps = event.angleDelta().y() // 120
		factor = 2**(-steps)
		self.move_plot("wheel", factor)

	def shortcuts(self):
		shortcuts_dict = {
			"w": lambda: self.move_plot("in"),
			"s": lambda: self.move_plot("out"),
			"a": lambda: self.move_plot("left"),
			"d": lambda: self.move_plot("right"),

			"Shift+w": lambda: self.move_plot("sin"),
			"Shift+s": lambda: self.move_plot("sout"),
			"Shift+a": lambda: self.move_plot("sleft"),
			"Shift+d": lambda: self.move_plot("sright"),
			
			# We cannot use shortcuts that are also used in the menu bar here, as these are global shortcuts on MacOS, as the MenuBar is always active
			"Alt+w": lambda: self.alter_plot(0, -1),
			"Alt+s": lambda: self.alter_plot(0, +1),
			"Alt+a": lambda: self.alter_plot(1, -1),
			"Alt+d": lambda: self.alter_plot(1, +1),
		}


		for keys, function in shortcuts_dict.items():
			QShortcut(keys, self.parent).activated.connect(function)

	@synchronized_d(locks["axs"])
	def on_range(self, xmin, xmax, index):
		axrange = self.axs["ax"][index].get_xlim()
		if xmax == xmin or xmax > axrange[1] or xmin < axrange[0]:
			return
		shift = (QApplication.keyboardModifiers() == Qt.KeyboardModifier.ShiftModifier)
		if self.onrange == "Assign" and (shift or main.config["fit_alwaysfit"]):
			xmiddle, xuncert = main.plotwidget.fit_peak(xmin, xmax, index)
			xpre = main.plotwidget.self.cache_positions[index]

			dict_ = {
				"x":		xmiddle,
				"error":	xuncert,
				"xpre":		xpre
			}

			reference = reference = main.plotwidget.get_series_reference(index[1])
			if reference["method"] == "Transition":
				qns = main.plotwidget.get_qns(reference["transition"])[index[0]]
				for i, qnu, qnl in zip(range(6), qns[0], qns[1]):
					dict_[f"qnu{i+1}"] = qnu
					dict_[f"qnl{i+1}"] = qnl
			elif reference["method"] == "List":
				qns = self.get_qns_list(reference, index=index)
				if qns is not None:
					for i, qnu, qnl in zip(range(6), qns[0], qns[1]):
						dict_[f"qnu{i+1}"] = qnu
						dict_[f"qnl{i+1}"] = qnl

			if not main.config["fit_alwaysassign"]:
				main.plotwidget.lastassignment = index, dict_
			elif main.plotwidget.check_blends(index, dict_):
				return
			else:
				main.plotwidget.assign(index, dict_)
		else:
			self.center = (xmin+xmax)/2
			self.width = xmax-xmin
			self.update_plot()

class TrendPlot(QWidget):
	def __init__(self, parent=None, get_data=None, get_init_level=None):
		super().__init__(parent)
		
		self.get_data = get_data
		self.get_init_level = get_init_level

		self.fit_data = pd.DataFrame(columns=pyckett.egy_dtypes.keys()).astype(pyckett.egy_dtypes)
		self.predicted_levels = None
		self.last_command = None
		self.red_egy_formula = "egy"
		self.last_level = None
		
		self.unique_qn = "qn1"
		self.init_gui()
		
	def init_gui(self):
		layout = QVBoxLayout()
		self.setLayout(layout)

		self.fig = figure.Figure(dpi=main.config["plot_dpi"])
		self.plot_canvas = FigureCanvas(self.fig)
		self.plot_canvas.setMinimumHeight(200)
		layout.addWidget(self.plot_canvas, 1)
		
		self.mpltoolbar = NavigationToolbar2QT(self.plot_canvas, self)
		layout.addWidget(self.mpltoolbar)
		
		self.axs = self.fig.subplots(2, sharex=True, gridspec_kw={"hspace": 0, "height_ratios": [5, 1]})
		self.axs[0].get_xaxis().set_visible(False)
		self.axs[1].get_yaxis().set_visible(False)
		self.axs[1].set_ylim(-1, 1)
		self.axs[1].get_xaxis().set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
		
		self.axs[0].callbacks.connect('xlim_changed', self.on_lims_change)
		self.axs[0].callbacks.connect('ylim_changed', self.on_lims_change)
		
		self.fit_line = matplotlib.collections.LineCollection(np.zeros(shape=(0,2,2)), colors=main.config["energylevelstrendwindow_fitcolor"])
		self.axs[0].add_collection(self.fit_line)
		
		self.fig.canvas.mpl_connect("motion_notify_event", self.on_hover)
		self.fig.canvas.mpl_connect("button_press_event", lambda x: self.on_hover(x, click=True))
		self.current_egy_data = None
		self.annot = self.axs[0].annotate("", xy=(0,0), xytext=(5,5), textcoords="offset points", color="black", ha="center", va="bottom", bbox=dict(boxstyle="round", fc="w"))
		self.annot.set_visible(False)
		
		self.egy_scatter = self.axs[0].scatter([], [], color=[], marker=".")
		self.ass_scatter = self.axs[0].scatter([], [], color=[], marker=".", zorder=10)
		self.res_scatter = self.axs[1].scatter([], [], color=[], marker=".")

		self.command_prompt = QQ(QLineEdit, placeholder="Enter commands (add, run, next, del, reset ...)")
		self.command_prompt.returnPressed.connect(self.on_command)
		layout.addWidget(self.command_prompt)
		
		self.best_levels_output = QQ(QTextEdit, readonly=True, maxHeight=160)
		self.best_levels_output.setFontFamily("Courier")
		layout.addWidget(self.best_levels_output)
		
		self.used_levels_output = QQ(QTextEdit, readonly=True, maxHeight=160)
		self.used_levels_output.setFontFamily("Courier")
		layout.addWidget(self.used_levels_output)

		for i in range(10):
			tmp = QShortcut(f"Ctrl+{i}", self.command_prompt)
			tmp.setContext(Qt.ShortcutContext.WidgetWithChildrenShortcut)
			tmp.activated.connect(lambda i=i: self.quick_assign(i))

	def on_lims_change(self, ax):
		xrange, yrange = self.axs[0].get_xlim(), self.axs[0].get_ylim()

		data = self.get_data().copy()
		data["red_egy"] = data.eval(self.red_egy_formula)
		query = f"{self.unique_qn} < {xrange[1]} and {self.unique_qn} > {xrange[0]} and red_egy < {yrange[1]} and red_egy > {yrange[0]}"
		data = data.query(query)

		offsets = data[[self.unique_qn, "red_egy"]].values
		color = main.config["energylevelswindow_defaultcolor"]
		self.egy_scatter.set(offsets=offsets, color=color)
		self.current_egy_data = data

	def quick_assign(self, i):
		if i == 0:
			i = 10

		if self.predicted_levels is None:
			main.notification("There are currently no predicted levels to be assigned.")
			return
		
		if i <= 0 or i > len(self.predicted_levels):
			main.notification("Specified index is out of bounds.")
			return
		
		self.add_energy_level(self.predicted_levels.iloc[i-1])

	def update_gui(self):
		columns_best = ["dist", "egy"] + [f"qn{i+1}" for i in range(main.config["series_qns"])]
		if self.predicted_levels is None:
			self.best_levels_output.setText("")
		else:
			self.best_levels_output.setText(f"Current best levels are:\n{self.predicted_levels[columns_best]}")
		
		columns = ["egy"] + [f"qn{i+1}" for i in range(main.config["series_qns"])]
		self.used_levels_output.setText(f"Currently used levels are:\n{self.fit_data.sort_values(self.unique_qn)[columns]}")

	def on_command(self, command=None):
		if command is None:
			command = self.command_prompt.text().lower().split()
		
		if not len(command):
			return
			
		elif command[0] == "del" and len(command) > 1:
			for qn in command[1:]:
				self.fit_data = self.fit_data.drop(self.fit_data[self.fit_data[self.unique_qn] == int(qn)].index)
			self.run_fit()
		
		elif command[0] == "reset":
			self.fit_data = self.fit_data.drop(self.fit_data.index)
			self.predicted_levels = None
			self.run_fit()
		
		elif command[0] == "add":
			if len(command) == 1:
				if self.predicted_levels is None:
					main.notification("Please specify a level to add, currently there are no predicted level.")
					return
				level_to_add = self.predicted_levels.iloc[0]
			else:
				qns = command[1:]
				query_string = " and ".join([f"qn{i+1} == {qn}" for i, qn in enumerate(qns)])

				level_to_add = self.get_data().query(query_string)
				if len(level_to_add) == 0:
					main.notification("Could not fing the specified level.")
					return
				
				if len(level_to_add) > 1:
					main.notification(f"Specified level was ambiguous, found {len(level_to_add)} matching levels.")
					return
				
				level_to_add = level_to_add.iloc[0]

			self.add_energy_level(level_to_add)
			
		elif command[0] == "next":
			if len(command) != 2:
				main.notification("Please specify exactly one argument to next.")
				return
			
			iterations = int(command[1])
			for i in range(iterations):
				response = self.add_energy_level(self.predicted_levels.iloc[0], dry=True)
				if response < 0:
					break
				
			self.run_fit()
		
		elif command[0] == "run":
			if len(command) == 1:
				self.run_fit()
			elif len(command) == 2:
				self.run_fit(next_level=int(command[1]))
			else:
				main.notification("Please specify exactly one or none arguments to run.")
				return
	
		elif len(command) == 1 and command[0] == ".":
			if self.last_command:
				self.on_command(command=self.last_command)
				return
			else:
				main.notification("No previous command specified")
				return
		
		else:
			main.notification("Command was not understood.")
			return
		
		self.last_command = command
		self.command_prompt.setText("")
	
	def run_fit(self, *args, **kwargs):
		try:
			resp = self.run_fit_core(*args, **kwargs)
		except Exception as E:
			main.notification(f"<span style='color:#eda711;'>WARNING</span>: The fit failed for the following reason:\n{E}")
			raise E
		if kwargs.get("dry"):
			return(resp)
		self.update_gui()
		self.plot_canvas.draw()
		return(resp)
		
	def run_fit_core(self, dry=False, next_level=None):
		return_value = 0
		red_egy_formula = self.red_egy_formula = main.config["energylevelstrendwindow_reducedenergy"] or "egy"
		
		assigned_xs = self.fit_data[self.unique_qn].values
		assigned_ys = self.fit_data.eval(red_egy_formula).values

		data = self.get_data().copy()
		data["red_egy"] = data.eval(red_egy_formula)

		prediction_proximity = main.config["energylevelstrendwindow_predictionproximity"]

		if not len(assigned_xs):
			main.notification("There is no suitable data for the fit.")
			self.predicted_levels = None
			
			qns = self.get_init_level()[:main.config["series_qns"]]
			query = " and ".join([f"qn{i+1} == {qn}" for i, qn in enumerate(qns)])
			init_levels = data.query(query).copy()
			
			init_levels["dist"] = 0
			init_levels["dist_abs"] = init_levels["dist"].abs()
			init_levels = init_levels.sort_values("dist_abs")
			self.predicted_levels = init_levels
			
			xinit = qns[int(self.unique_qn[2:])-1]
			plot_xrange = (xinit - prediction_proximity, xinit + prediction_proximity)
			
			if len(init_levels):
				ymean = init_levels["red_egy"].mean()
				plot_yrange = (ymean - 1, ymean + 1)
			else:
				plot_yrange = (-1, 1)
			
			empty_offsets = np.empty((0, 2))
			self.ass_scatter.set(offsets=empty_offsets)
			self.res_scatter.set(offsets=empty_offsets)
			
			self.axs[0].set_xlim(*plot_xrange)
			self.axs[0].set_ylim(*plot_yrange)
			return -1
		
		xmin, xmax = assigned_xs.min(), assigned_xs.max()
		
		if next_level is not None:
			next_level = next_level
		elif self.last_level is not None:
			next_level = self.last_level + 1
		else:
			next_level = xmax + 1
		
		fit_proximity = main.config["energylevelstrendwindow_fitproximity"]
		if fit_proximity:
			mask = np.abs(assigned_xs - next_level) <= fit_proximity
			xs_used_in_fit = assigned_xs[mask]
			ys_used_in_fit = assigned_ys[mask]
		else:
			xs_used_in_fit = assigned_xs
			ys_used_in_fit = assigned_ys
		
		if not len(xs_used_in_fit):
			main.notification("There is no suitable data for the fit.")
			self.predicted_levels = None
			return -1
		
		rank = min( max(0, len(xs_used_in_fit) - 2), main.config["energylevelstrendwindow_maxrank"])
		polynom_from_fit = np.polyfit(xs_used_in_fit, ys_used_in_fit, rank)
		
		predicted_xs = np.arange(next_level - fit_proximity, next_level + fit_proximity + 1, 0.2)
		predicted_ys = np.polyval(polynom_from_fit, predicted_xs)

		next_energy = np.polyval(polynom_from_fit, next_level)
		
		max_deviation = main.config["energylevelstrendwindow_maxdeviation"]
		possible_levels = data.query(f"{self.unique_qn} == {next_level} and abs(red_egy - {next_energy}) < {max_deviation}")
		
		if not len(possible_levels):
			main.notification("Did not find any possible next levels.")
			return_value = -1
			self.predicted_levels = None
		
		else:
			possible_levels = possible_levels.copy()
			possible_levels["dist"] = possible_levels["red_egy"] - next_energy
			possible_levels["dist_abs"] = possible_levels["dist"].abs()
			possible_levels = possible_levels.sort_values("dist_abs")
			self.predicted_levels = possible_levels.head(10)
		
		if dry:
			return return_value
		
		segs = np.array(((predicted_xs[:-1], predicted_xs[1:]), (predicted_ys[:-1], predicted_ys[1:]))).T
		color = main.config["energylevelstrendwindow_fitcolor"]
		self.fit_line.set(segments=segs, color=color)
		
		offsets = np.vstack((assigned_xs, assigned_ys)).T
		self.ass_scatter.set(offsets=offsets, color=color)
		
		offsets[:, 1] = 0
		color = "green"
		self.res_scatter.set(offsets=offsets, color=color)
		
		if main.config["energylevelstrendwindow_plotshowall"]:
			plot_xrange = (xmin - prediction_proximity, xmax + prediction_proximity)
		else:
			plot_xrange = (next_level - prediction_proximity, next_level + prediction_proximity)
		
		mask = (predicted_xs >= plot_xrange[0]) & (predicted_xs <= plot_xrange[1])
		plot_yrange = (predicted_ys[mask].min(), predicted_ys[mask].max())
		if plot_yrange[0] == plot_yrange[1]:
			plot_yrange = (plot_yrange[0] - 1, plot_yrange[0] + 1)
		plot_yrange_dist = plot_yrange[1] - plot_yrange[0]
		margin = main.config["plot_ymargin"]
		plot_yrange = (plot_yrange[0] - plot_yrange_dist * margin, plot_yrange[1] + plot_yrange_dist * margin)
		
		self.axs[0].set_xlim(*plot_xrange)
		self.axs[0].set_ylim(*plot_yrange)
		
		self.axs[1].set_ylim(-1, 1)
		
		# The total fit data is done in on_lims_change -> Allows to pan and zoom interactively
		
		return return_value
	
	def add_energy_level(self, energy_level_row, dry=False):
		new_qn = energy_level_row[self.unique_qn]
		self.last_level = new_qn
		self.fit_data = self.fit_data.drop(self.fit_data[self.fit_data[self.unique_qn] == new_qn].index)
		self.fit_data = self.fit_data.append(energy_level_row)
		resp = self.run_fit(dry=dry)
		return(resp)

	def on_hover(self, event, click=False):
		if event.inaxes == self.axs[0] and self.current_egy_data is not None:
			cont, ind = self.egy_scatter.contains(event)
			if cont:
				indices = ind["ind"]
				noq = main.config["series_qns"]
				tmp_levels = self.current_egy_data.iloc[indices]
				
				text = []
				for i, row in tmp_levels.iterrows():
					text.append(",".join(str(int(row[f"qn{i+1}"])) for i in range(noq)))
				text = "\n".join(text)
				
				self.annot.xy = self.egy_scatter.get_offsets()[indices[0]]
				self.annot.set(text=text, visible=True)
				if click:
					level_to_add = tmp_levels.iloc[0]
					self.add_energy_level(level_to_add)
			else:
				self.annot.set_visible(False)
			self.fig.canvas.draw()
			

class NotificationsBox(QWidget):
	def __init__(self):
		super().__init__()
		self.bg_color = QColor("#a5aab3")
		self.messages = []
		self.setWindowFlags(
			Qt.WindowType.Window | Qt.WindowType.Tool | Qt.WindowType.FramelessWindowHint |
			Qt.WindowType.WindowStaysOnTopHint | Qt.WindowType.X11BypassWindowManagerHint)

		self.setAttribute(Qt.WidgetAttribute.WA_NoSystemBackground, True)
		self.setAttribute(Qt.WidgetAttribute.WA_TranslucentBackground, True)

		self.setMinimumHeight(80)
		self.setMinimumWidth(300)
		self.setMaximumWidth(300)

		self.layout = QVBoxLayout()
		self.setLayout(self.layout)

		self.setStyleSheet("""
			color: white;
			background-color: #bf29292a;
		""")

		self._desktop = QApplication.instance().primaryScreen()
		startPos = QPoint(self._desktop.geometry().width() - self.width() - 10, 10)
		self.move(startPos)

	def paintEvent(self, event=None):
		painter = QPainter(self)

		painter.setOpacity(0.5)
		painter.setPen(QPen(self.bg_color))
		painter.setBrush(self.bg_color)
		painter.drawRect(self.rect())

	def add_message(self, text):
		label = QLabel(text)
		label.setWordWrap(True)
		label.setStyleSheet("""
			padding: 5px;
		""")

		self.layout.addWidget(label)
		self.messages.append(label)
		self.timer = QTimer(self)
		self.timer.setSingleShot(True)
		self.timer.timeout.connect(self.unshow_message)
		self.timer.start(main.config["flag_notificationtime"])

		self.show()

		self.timer2 = QTimer(self)
		self.timer2.setSingleShot(True)
		self.timer2.timeout.connect(self.adjustSize)
		self.timer2.start(0)


	def unshow_message(self):
		label = self.messages.pop()
		label.hide()
		label.deleteLater()
		if not self.messages:
			self.hide()
		self.adjustSize()

class CustomError(Exception):
	pass

class CDict(dict):
	def __init__(self, callback, *args, **kwargs):
		self.callback = callback
		return(super().__init__(*args, **kwargs))

	def __setitem__(self, key, value):
		if isinstance(value, dict):
			value = CDict(self.callback, value)
		super().__setitem__(key, value)
		self.callback()

	def update(self, items):
		if isinstance(items, dict):
			items = items.items()
		for key, value in items:
			self[key] = value

class Config(dict):
	def __init__(self, signal):
		super().__init__()
		self.signal = signal
		self.signal.connect(self.callback)
		self.callbacks = pd.DataFrame(columns=["id", "key", "widget", "function"], dtype="object").astype({"id": np.uint})


	def __setitem__(self, key, value, widget=None):
		super().__setitem__(key, value)
		self.signal.emit((key, value, widget))

	def callback(self, args):
		key, value, widget = args
		if widget:
			callbacks_widget = self.callbacks.query(f"key == @key and widget != @widget")
		else:
			callbacks_widget = self.callbacks.query(f"key == @key")
		for i, row in callbacks_widget.iterrows():
			row["function"]()

	def register(self, keys, function):
		if not isinstance(keys, (tuple, list)):
			keys = [keys]
		for key in keys:
			id = 0
			df = self.callbacks
			df.loc[len(df), ["id", "key", "function"]] = id, key, function

	def register_widget(self, key, widget, function):
		ids = set(self.callbacks["id"])
		id = 1
		while id in ids:
			id += 1
		df = self.callbacks
		df.loc[len(df), ["id", "key", "function", "widget"]] = id, key, function, widget
		widget.destroyed.connect(lambda x, id=id: self.unregister_widget(id))

	def unregister_widget(self, id):
		self.callbacks.drop(self.callbacks[self.callbacks["id"] == id].index, inplace=True)

def is_dark_theme():
	return(QApplication.styleHints().colorScheme() == Qt.ColorScheme.Dark)

class Color(str):
	def __new__(cls, color):
		cls.validate_color(cls, color)
		return super().__new__(cls, color)

	def __assign__(self, color):
		self.validate_color(color)
		return super().__new__(color)

	def validate_color(self, color):
		match = re.search(r'^#(?:[0-9a-fA-F]{3}?){1,2}$|^#(?:[0-9a-fA-F]{8}?)$', color)
		if match:
			if len(color) == 9 and color[-2:] == "ff":
				color = color[:-2]
			return(color)
		else:
			raise CustomError(f"Invalid Color: '{color}' is not a valid color.")

class SignalClass(QObject):
	fileschanged    = pyqtSignal()
	assignment      = pyqtSignal()
	updateplot      = pyqtSignal()
	drawplot        = pyqtSignal()
	createdplots    = pyqtSignal()
	blwfit          = pyqtSignal()
	peakfinderstart = pyqtSignal()
	peakfinderend   = pyqtSignal()
	calibrationstart= pyqtSignal()
	calibrationend  = pyqtSignal()
	overlapend      = pyqtSignal()
	overlapindicator= pyqtSignal(str)
	fitindicator    = pyqtSignal(str)
	setindicator    = pyqtSignal(str)
	writelog        = pyqtSignal(str)
	writehover      = pyqtSignal(str)
	notification    = pyqtSignal(str)
	updateconfig    = pyqtSignal(tuple)
	def __init__(self):
		super().__init__()

class QSpinBox(QSpinBox):
	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		# AdaptiveDecimalStepType is not implemented in earlier versions of PyQt5
		try:
			self.setStepType(QAbstractSpinBox.StepType.AdaptiveDecimalStepType)
		except:
			pass

	def setSingleStep(self, value):
		self.setStepType(QAbstractSpinBox.StepType.DefaultStepType)
		super().setSingleStep(value)

	def setValue(self, value):
		if value < -2147483647 or value > 2147483647:
			value = 0
		return super().setValue(value)

	def setRange(self, min, max):
		min = min if not min is None else -2147483647
		max = max if not max is None else +2147483647
		return super().setRange(min, max)

class QDoubleSpinBox(QDoubleSpinBox):
	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		self.setDecimals(20)
		# AdaptiveDecimalStepType is not implemented in earlier versions of PyQt5
		try:
			self.setStepType(QAbstractSpinBox.StepType.AdaptiveDecimalStepType)
		except:
			pass

	def setSingleStep(self, value):
		self.setStepType(QAbstractSpinBox.StepType.DefaultStepType)
		super().setSingleStep(value)

	def textFromValue(self, value):
		if value and abs(np.log10(abs(value))) > 5:
			return(f"{value:.2e}")
		else:
			return(f"{value:.10f}".rstrip("0").rstrip("."))

	def valueFromText(self, text):
		return(np.float64(text))

	def setRange(self, min, max):
		min = min if not min is None else -np.inf
		max = max if not max is None else +np.inf
		return super().setRange(min, max)

	def validate(self, text, position):
		try:
			np.float64(text)
			return(QValidator.State(2), text, position)
		except ValueError:
			if text.strip() in ["+", "-", ""]:
				return(QValidator.State(1), text, position)
			elif re.match(r"^[+-]?\d+\.?\d*[Ee][+-]?\d?$", text):
				return(QValidator.State(1), text, position)
			else:
				return(QValidator.State(0), text, position)

	def fixup(self, text):
		tmp = re.search(r"[+-]?\d+\.?\d*", text)
		if tmp:
			return(tmp[0])
		else:
			return(str(0))

class QTableWidget(QTableWidget):
	def keyPressEvent(self, event):
		if not csv_copypaste(self, event):
			super().keyPressEvent(event)

class QTableView(QTableView):
	def keyPressEvent(self, event):
		if not csv_copypaste(self, event):
			super().keyPressEvent(event)

class CustomTableModel(QAbstractTableModel):
	def __init__(self, data, headers=None, hiddencolumns=[], editable=True):
		super().__init__()

		self.data = data
		self.headers = headers if headers else data.columns
		self.hiddencolumns = hiddencolumns
		self.editable = editable

	def data(self, index, role):
		if role == Qt.ItemDataRole.DisplayRole or role == Qt.ItemDataRole.EditRole:
			value = self.data.iloc[index.row(), index.column()]
			dtype = self.data[self.data.columns[index.column()]].dtypes

			if isinstance(value, str):
				return(value)
			elif isinstance(value, (np.integer, int)):
				if value == pyckett.SENTINEL:
					return("")
				else:
					if role == Qt.ItemDataRole.EditRole:
						return(str(value))
					else:
						return(f"{{:{main.config['flag_tableformatint']}}}".format(value))
			elif np.isnan(value):
				return("")
			else:
				if role == Qt.ItemDataRole.EditRole:
					return(str(value))
				else:
					return(f"{{:{main.config['flag_tableformatfloat']}}}".format(value))

	def rowCount(self, index):
		return(self.data.shape[0])

	def columnCount(self, index):
		return(self.data.shape[1]-len(self.hiddencolumns))

	def headerData(self, section, orientation, role):
		if role == Qt.ItemDataRole.DisplayRole:
			if orientation == Qt.Orientation.Horizontal:
				return str(self.headers[section])

			if orientation == Qt.Orientation.Vertical:
				if section >= len(self.data.index):
					return ""
				return str(self.data.index[section])

	def flags(self, index):
		if self.editable:
			return(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable | Qt.ItemFlag.ItemIsEditable)
		else:
			return(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable)

	def update(self):
		self.layoutChanged.emit()

	def setData(self, index, value, role):
		if not index.isValid():
			return False
		if role != Qt.ItemDataRole.EditRole:
			return False

		row = index.row()
		if row < 0 or row >= len(self.data):
			return False
		column = index.column()
		if column < 0 or column >= self.data.columns.size:
			return False

		dtype = self.data[self.data.columns[index.column()]].dtypes
		if np.issubdtype(dtype, np.number):
			try:
				if np.issubdtype(dtype, np.integer):
					value = np.int64(value)
				else:
					value = np.float64(value)
			except ValueError:
				if np.issubdtype(dtype, np.integer):
					value = pyckett.SENTINEL
				else:
					value = np.nan

		self.data.iloc[row, column] = value
		self.dataChanged.emit(index, index)
		return True




##
## Global Functions
##
def QQ(widgetclass, config_key=None, **kwargs):
	widget = widgetclass()

	if "range" in kwargs:
		widget.setRange(*kwargs["range"])
	if "maxWidth" in kwargs:
		widget.setMaximumWidth(kwargs["maxWidth"])
	if "maxHeight" in kwargs:
		widget.setMaximumHeight(kwargs["maxHeight"])
	if "minWidth" in kwargs:
		widget.setMinimumWidth(kwargs["minWidth"])
	if "minHeight" in kwargs:
		widget.setMinimumHeight(kwargs["minHeight"])
	if "color" in kwargs:
		widget.setColor(kwargs["color"])
	if "text" in kwargs:
		widget.setText(kwargs["text"])
	if "options" in kwargs:
		options = kwargs["options"]
		if isinstance(options, dict):
			for key, value in options.items():
				widget.addItem(key, value)
		else:
			for option in kwargs["options"]:
				widget.addItem(option)
	if "width" in kwargs:
		widget.setFixedWidth(kwargs["width"])
	if "height" in kwargs:
		widget.setFixedHeight(kwargs["height"])
	if "tooltip" in kwargs:
		widget.setToolTip(kwargs["tooltip"])
	if "placeholder" in kwargs:
		widget.setPlaceholderText(kwargs["placeholder"])
	if "singlestep" in kwargs:
		widget.setSingleStep(kwargs["singlestep"])
	if "wordwrap" in kwargs:
		widget.setWordWrap(kwargs["wordwrap"])
	if "align" in kwargs:
		widget.setAlignment(kwargs["align"])
	if "rowCount" in kwargs:
		widget.setRowCount(kwargs["rowCount"])
	if "columnCount" in kwargs:
		widget.setColumnCount(kwargs["columnCount"])
	if "move" in kwargs:
		widget.move(*kwargs["move"])
	if "default" in kwargs:
		widget.setDefault(kwargs["default"])
	if "textFormat" in kwargs:
		widget.setTextFormat(kwargs["textFormat"])
	if "checkable" in kwargs:
		widget.setCheckable(kwargs["checkable"])
	if "shortcut" in kwargs:
		widget.setShortcut(kwargs["shortcut"])
	if "parent" in kwargs:
		widget.setParent(kwargs["parent"])
	if "completer" in kwargs:
		widget.setCompleter(kwargs["completer"])
	if "hidden" in kwargs:
		widget.setHidden(kwargs["hidden"])
	if "visible" in kwargs:
		widget.setVisible(kwargs["visible"])
	if "stylesheet" in kwargs:
		widget.setStyleSheet(kwargs["stylesheet"])
	if "enabled" in kwargs:
		widget.setEnabled(kwargs["enabled"])
	if "items" in kwargs:
		for item in kwargs["items"]:
			widget.addItem(item)
	if "readonly" in kwargs:
		widget.setReadOnly(kwargs["readonly"])
	if "prefix" in kwargs:
		widget.setPrefix(kwargs["prefix"])

	if widgetclass in [QSpinBox, QDoubleSpinBox]:
		setter = widget.setValue
		changer = widget.valueChanged.connect
		getter = widget.value
	elif widgetclass == QCheckBox:
		setter = widget.setChecked
		changer = widget.stateChanged.connect
		getter = widget.isChecked
	elif widgetclass == QPlainTextEdit:
		setter = widget.setPlainText
		changer = widget.textChanged.connect
		getter = widget.toPlainText
	elif widgetclass == QLineEdit:
		setter = widget.setText
		changer = widget.textChanged.connect
		getter = widget.text
	elif widgetclass == QAction:
		setter = widget.setChecked
		changer = widget.triggered.connect
		getter = widget.isChecked
	elif widgetclass == QPushButton:
		setter = widget.setDefault
		changer = widget.clicked.connect
		getter = widget.isDefault
	elif widgetclass == QToolButton:
		setter = widget.setChecked
		changer = widget.clicked.connect
		getter = widget.isChecked
	elif widgetclass == QComboBox:
		setter = widget.setCurrentText
		changer = widget.currentTextChanged.connect
		getter = widget.currentText
	else:
		return widget

	if "value" in kwargs:
		setter(kwargs["value"])
	if config_key:
		setter(main.config[config_key])
		changer(lambda x=None, key=config_key: main.config.__setitem__(key, getter(), widget))
		main.config.register_widget(config_key, widget, lambda: setter(main.config[config_key]))
	if "change" in kwargs:
		changer(kwargs["change"])
	if "changes" in kwargs:
		for change in kwargs["changes"]:
			changer(change)

	return widget

def bin_data(dataframe, binwidth, range):
	length = len(dataframe)
	
	dataframe.loc[:,"bin"] = (dataframe.loc[:,"x"]-range[0]) // binwidth

	# For assignments (lin_df) as they do not have an intensity
	if "y" not in dataframe:
		dataframe = dataframe.loc[dataframe.drop_duplicates(("bin", "filename"), keep="last").sort_values(["x"]).index]
	else:
		dataframe = dataframe.loc[dataframe.sort_values("y").drop_duplicates(("bin", "filename"), keep="last").sort_values(["x"]).index]
	return(dataframe)

def csv_copypaste(self, event):
	if event.key() == Qt.Key.Key_C and (event.modifiers() == Qt.KeyboardModifier.ControlModifier):
		cells = sorted(self.selectedIndexes())
		output = []
		i = 0

		while i < len(cells):
			tmp = []
			row = cells[i].row()
			while i < len(cells) and cells[i].row() == row:
				tmp.append(cells[i].data())
				i += 1
			output.append("\t".join(map(str, tmp)))
		output = "\n".join(output)
		QApplication.clipboard().setText(output)

	elif event.key() == Qt.Key.Key_V and (event.modifiers() == Qt.KeyboardModifier.ControlModifier):
		if QAbstractItemView.EditTrigger.NoEditTriggers == self.editTriggers():
			return
		cells = sorted(self.selectedIndexes())
		if not cells:
			return
		text = QApplication.clipboard().text()
		data = [row.split("\t") for row in text.split("\n")]
		i_0, j_0 = cells[0].row(), cells[0].column()

		for i, row in enumerate(data):
			j_hidden = 0
			for j, value in enumerate(row):
				while self.isColumnHidden(j_0+j+j_hidden):
					j_hidden += 1
				self.model().setData(self.model().index(i_0+i, j_0+j+j_hidden), value, Qt.ItemDataRole.EditRole)
	else:
		return False
	return True

def trgb_to_rgbt(color):
	if len(color) == 9:
		color = f"#{color[3:]}{color[1:3]}"
	return(color)

def rgbt_to_trgb(color):
	if len(color) == 9:
		color = f"#{color[-2:]}{color[1:-2]}"
	return(color)

def encodebytes(bytes):
	return ", ".join([str(x) for x in bytes])

def decodebytes(code):
	return bytearray([int(x) for x in code.split(", ")])

def create_colors(dataframe, files={}, xpos=None, lin=False):
	if lin:
		files = files.copy()
		files["__lin__"] = {"color": main.config["color_lin"]}
	if len(files)==1 and xpos==None:
		return(files[list(files.keys())[0]].get("color", "#ffffff"))
	if not isinstance(dataframe, pd.DataFrame):
		return([])

	dataframe.reset_index(drop=True, inplace=True)
	tmp_colors = {file: files[file].get("color", "#ffffff") for file in files.keys()}
	filenames = dataframe["filename"]
	colors = filenames.replace(tmp_colors).copy()

	if xpos:
		colors[dataframe["x"] == xpos] = main.config["color_cur"]

	return(colors)

def column_to_numeric(val, force_int=False):
	val = val.strip()
	if val == "" or val == ":" or val == len(val)*"*":
		if force_int:
			return(pyckett.SENTINEL)
		else:
			return(np.nan)
	elif val[0].isalpha():
		val = str(ord(val[0].upper())-55)+val[1:]

	if force_int:
		return(np.int64(val))
	else:
		return(np.float64(val))

def symmetric_ticklabels(ticks):
	tick_labels = []
	for a, o in zip(ticks, ticks[::-1]):
		if not (np.isfinite(a) and np.isfinite(o)):
			continue
		dec_a = len(f"{a:{main.config['flag_xformatfloat']}}".rstrip("0").split(".")[1])
		dec_o = len(f"{o:{main.config['flag_xformatfloat']}}".rstrip("0").split(".")[1])
		if dec_a == dec_o:
			tick_labels.append(f"{a:{main.config['flag_xformatfloat']}}".rstrip("0").rstrip("."))
		else:
			trailing_zeros = 4 - max(dec_a, dec_o)
			tick = f"{a:{main.config['flag_xformatfloat']}}"[:-trailing_zeros] if trailing_zeros else f"{a:{main.config['flag_xformatfloat']}}"
			tick_labels.append(tick)
	return(tick_labels)

def except_hook(cls, exception, traceback):
	if issubclass(cls, KeyboardInterrupt):
		sys.exit(0)

	sys.__excepthook__(cls, exception, traceback)
	with open(llwpfile(".err"), "a+", encoding="utf-8") as file:
		time_str = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
		file.write(f"{time_str}: \n{exception}\n{''.join(tb.format_tb(traceback))}\n\n")
	try:
		if not main.new_df.empty:
			main.save_lines_lin(llwpfile(".lin"), force_noappend=True, force_lin=True)
		if main.config["flag_debug"]:
			main.self.notification(f"{exception}\n{''.join(tb.format_tb(traceback))}")
	except Exception as E:
		pass

def earlyreturn(ownid, lastid):
	if ownid != lastid:
		raise CustomError()

def commandline(showdialog=True):
	if showdialog:
		dialog = ConsoleDialog()
		dialog.exec()
		if dialog.result() != 1:
			return

	title, command = main.config["commandlinedialog_commands"][main.config["commandlinedialog_current"]]
	if not command.strip():
		return
	message = []
	old_stdout = sys.stdout
	red_output = sys.stdout = io.StringIO()
	try:
		exec(command)
	except Exception as E:
		message.append(f"<span style='color:#eda711;'>WARNING</span>: Executing the code raised an error: {str(E)}")
		raise
	finally:
		sys.stdout = old_stdout

	message.append("\n".join([f">>> {line}" for line in command.split("\n")]))
	message.append(red_output.getvalue())
	main.notification("\n".join(message))

def send_mail_to_author():
	webbrowser.open(f"mailto:bonah@ph1.uni-koeln.de?subject={APP_TAG}")

def restart():
	fname = llwpfile(".llwp")
	main.saveproject(fname)
	main.saveoptions()
	if not main.new_df.empty:
		main.save_lines_lin(llwpfile(".lin"), force_noappend=True, force_lin=True)
	os.execl(sys.executable, sys.executable, __file__, fname)

def lineshape(shape, derivative, *args):
	if shape == "Gauss":
		x, x_0, amp, width = args
		width = width/(2*np.sqrt(2*np.log(2)))
		if width == 0:
			return [0 if i!=x_0 else np.inf for i in x]
		ys = np.exp(-(x-x_0)**2/(2*width**2))/(width*(2*np.pi)**0.5)
	elif shape == "Lorentz":
		x, x_0, amp, width = args
		width = width/2
		if width == 0:
			return [0 if i!=x_0 else np.inf for i in x]
		ys = 1/(np.pi*width*(1+((x-x_0)/width)**2))
	else: #Voigt
		x, x_0, amp, gauss, lorentz = args
		gauss = gauss/(2*np.sqrt(2*np.log(2)))
		lorentz = lorentz/2
		if gauss == 0 and lorentz == 0:
			return [0 if i!=x_0 else np.inf for i in x]
		ys = special.voigt_profile(x-x_0, gauss, lorentz)

	for i in range(0, int(derivative)):
		ys = np.gradient(ys, edge_order=2)
	if derivative%2 == 0 and derivative != 0 and derivative%4 != 0:
		ys = -ys
	ymax = np.max(ys) if np.isfinite(ys).any() else 1
	if not np.isfinite(ymax) or ymax == 0:
		ymax = 1
	ys = amp*ys/ymax
	return(ys)

def addemptyrow_inplace(df, model=None):
	df.reset_index(drop=True, inplace=True)
	dtypes = df.dtypes
	newvalues = []
	for dtype in dtypes:
		if dtype == np.float64:
			newvalues.append(np.nan)
		elif dtype == np.int64:
			newvalues.append(pyckett.SENTINEL)
		else:
			newvalues.append("")
	df.loc[len(df.index)] = newvalues
	if model:
		model.update()

def exp_to_df(fname, sep="\t", xcolumn=0, ycolumn=1, sort=True):
	data = pd.read_csv(fname, sep=sep, dtype=np.float64, header=None, engine="c", comment="#")
	column_names = [i for i in range(len(data.columns))]

	column_names[xcolumn if xcolumn in column_names else 0] = "x"
	column_names[ycolumn if ycolumn in column_names else 1] = "y"

	data.columns = column_names
	data = data[["x", "y",]]
	data = data.dropna()
	data["filename"] = fname

	return(data)

def shared_labels(fig, xlabel, ylabel, xlabelpad=15, ylabelpad=0, **kwargs):
	ax = fig.add_subplot(111, frameon=False)
	for side in ("top", "bottom", "left", "right"):
		ax.spines[side].set_color('none')
	ax.tick_params(axis="both", labelcolor='#fff0', top=False, bottom=False, left=False, right=False, zorder=-100)
	ax.set_xlabel(xlabel, labelpad=xlabelpad, **kwargs)
	ax.set_ylabel(ylabel, labelpad=ylabelpad, **kwargs)
	return(ax)

def llwpfile(extension):
	home = os.path.expanduser("~")
	llwpfolder = os.path.join(home, ".llwp")
	
	if not os.path.isdir(llwpfolder):
		os.mkdir(llwpfolder)

	return(os.path.join(llwpfolder, extension))


##
## Global Variables
##

exp_dtypes = {
  'x':			np.float64,
  'y':			np.float64,
}

cat_dtypes = pyckett.cat_dtypes
lin_dtypes = pyckett.lin_dtypes

config_specs = {
	# Format is: [default value, class]
	"layout_mpltoolbar":					[False, bool],

	"color_exp":							["#000000", Color],
	"color_lin":							["#ff38fc", Color],
	"color_cat":							["#d91e6f", Color],
	"color_cur":							["#71eb34", Color],
	"color_fit":							["#bc20e3", Color],

	"fit_error":							[0.05, float],
	"fit_function":							["Polynom", str],
	"fit_alwaysfit":						[True, bool],
	"fit_alwaysassign":						[True, bool],
	"fit_polynomrank":						[5, int],
	"fit_polynommaxrank":					[20, int],
	"fit_uncertaintystep":					[0.01, float],
	"fit_xpoints":							[1000, int],
	"fit_offset":							[True, bool],
	"fit_comment":							["", str],
	"fit_clipboard":						[True, bool],
	"fit_amplitudedirection":				[1, int],

	"plot_dpi":								[100, float],
	"plot_annotations_dict":				[{"x": 1, "y": 1, "horizontalalignment": "right", "verticalalignment": "top"}, dict],
	"plot_font_dict":						[{"size":10}, dict],
	"plot_matplotlibkwargs":				[{"hspace": 0, "wspace": 0}, dict],
	"plot_coupled":							[True, bool],
	"plot_ymargin":							[0.1, float],
	"plot_annotate":						[True, bool],
	"plot_hover_cutoff":					[20, float],
	"plot_rows":							[5, int],
	"plot_cols":							[1, int],
	"plot_ticks":							[3, int],
	"plot_offsetticks":						[1, int],
	"plot_yscale":							["Per Plot", str],
	"plot_expcat_factor":					[1, float],
	"plot_expcat_exponent":					[10, int],
	"plot_yscale_min":						[-100, float],
	"plot_yscale_max":						[300, float],
	"plot_offset":							[0, float],
	"plot_width":							[100, float],
	"plot_bins":							[4000, int],
	"plot_skipbinning":						[1000, int],
	"plot_expasstickspectrum":				[False, bool],
	"plot_relativegoto":					[False, bool],

	"series_qns":							[3, int],
	"series_annotate_xs":					[False, bool],
	"series_annotate_fmt":					[".4f", str],
	"series_blendwidth":					[1, float],
	"series_blendminrelratio":				[0, float],
	"series_blenddialog":					[True, bool],
	"series_currenttab":					[1, int],
	"series_references":					[[], list],

	"flag_automatic_draw":					[True, bool],
	"flag_appendonsave":					[True, bool],
	"flag_hidecatalog":						[False, bool],
	"flag_xcolumn":							[0, int],
	"flag_ycolumn":							[1, int],
	"flag_separator":						[9, int],
	"flag_debug":							[False, bool],
	"flag_alwaysshowlog":					[True, bool],
	"flag_extensions":						[{"exp": [".csv"], "cat": [".cat"], "lin": [".lin"], "project": [".llwp"]}, dict],
	"flag_autosetqns":						[True, bool],
	"flag_predictionformats":				[{}, dict],
	"flag_assignmentformats":				[{}, dict],
	"flag_assignmentsavefmt":				[{}, dict],
	"flag_loadfilesthreaded":				[True, bool],
	"flag_shownotification":				[True, bool],
	"flag_notificationtime":				[2000, int],
	"flag_showmainplotcontrols":			[True, bool],
	"flag_showmainplotwidth":				[True, bool],
	"flag_showmainplotrowscols":			[True, bool],
	"flag_showmainplotposition":			[True, bool],
	"flag_referencenumberlocked":			[True, bool],
	"flag_logmaxrows":						[10000, int],
	"flag_tableformatint":					[".0f", str],
	"flag_tableformatfloat":				[".2f", str],
	"flag_xformatfloat":					[".4f", str],
	"flag_autosave":						[120, int],
	"flag_protplotautowidth":				[True, bool],

	"pipe_current":							[0, int],
	"pipe_commands":						[[], list],

	"blendedlineswindow_lineshape":			["Gauss", str],
	"blendedlineswindow_derivative":		[0, int],
	"blendedlineswindow_transparency":		[0.2, float],
	"blendedlineswindow_maxfwhm":			[10, float],
	"blendedlineswindow_polynom":			[0, int],
	"blendedlineswindow_fixedwidth":		[False, bool],
	"blendedlineswindow_showbaseline":		[True, bool],
	"blendedlineswindow_xpoints":			[1000, int],
	"blendedlineswindow_color_total":		["#3d5dff", Color],
	"blendedlineswindow_color_points":		["#ff3352", Color],
	"blendedlineswindow_color_baseline":	["#f6fa14", Color],
	"blendedlineswindow_autopositionpeaks":	[True, bool],

	"lineshapewindow_lineshape":			["Gauss", str],
	"lineshapewindow_derivative":			[0, int],
	"lineshapewindow_scaling_factor":		[1, float],
	"lineshapewindow_scaling_exponent":		[10, int],
	"lineshapewindow_gauss":				[1, float],
	"lineshapewindow_lorentz":				[1, float],
	"lineshapewindow_xpoints":				[10000, int],

	"seriesfinderwindow_start":				["", str],
	"seriesfinderwindow_stop":				["", str],
	"seriesfinderwindow_results":			[10, int],
	"seriesfinderwindow_condition":			["", str],
	"seriesfinderwindow_atype":				[True, bool],
	"seriesfinderwindow_btype":				[True, bool],
	"seriesfinderwindow_ctype":				[True, bool],
	"seriesfinderwindow_onlyunassigned":	[True, bool],

	"seriesfitwindow_function":				["", str],
	"seriesfitwindow_maxprediction":		[100, int],
	"seriesfitwindow_series":				[{}, dict],
	"seriesfitwindow_greedy":				[True, bool],

	"peakfinderwindow_peakcolor":			["#4287f5", Color],
	"peakfinderwindow_kwargs":				[{}, dict],
	"peakfinderwindow_onlyunassigned":		[True, bool],
	"peakfinderwindow_width":				[1, float],
	"peakfinderwindow_maxentries":			[1000, int],

	"qnsdialog_width":						[1000, int],
	"qnsdialog_height":						[500, int],

	"commandlinedialog_width":				[500, int],
	"commandlinedialog_height":				[250, int],
	"commandlinedialog_commands":			[[], list],
	"commandlinedialog_current":			[1, int],

	"residualswindow_defaultcolor":			["#000000", Color],
	"residualswindow_query":				["", str],
	"residualswindow_colorinput":			["", str],
	"residualswindow_xvariable":			["", str],
	"residualswindow_yvariable":			["", str],
	"residualswindow_autoscale":			[True, bool],
	"residualswindow_blends":				[False, bool],

	"energylevelswindow_defaultcolor":		["#000000", Color],
	"energylevelswindow_query":				["", str],
	"energylevelswindow_colorinput":		["", str],
	"energylevelswindow_xvariable":			["", str],
	"energylevelswindow_yvariable":			["", str],
	"energylevelswindow_autoscale":			[True, bool],
	
	"energylevelstrendwindow_fitcolor":					["#ff1a74", Color],
	"energylevelstrendwindow_reducedenergy":		 	["", str],
	"energylevelstrendwindow_maxrank":					[4, int],
	"energylevelstrendwindow_maxdeviation":				[1, float],
	"energylevelstrendwindow_fitproximity":				[20, float],
	"energylevelstrendwindow_predictionproximity":		[5, float],
	"energylevelstrendwindow_egycatconversionfactor":	[29979.2458, float],
	"energylevelstrendwindow_egyquery":					["", str],
	"energylevelstrendwindow_plotshowall":				[False, bool],
	"energylevelstrendwindow_egyqns":					[True, bool],

	"reportwindow_blends":					[True, bool],
	"reportwindow_query":					["", str],

	"figurewindow_width":					[14, float],
	"figurewindow_height":					[8.65, float],
	"figurewindow_unit":					["cm", str],
	"figurewindow_dpi":						[300, float],
	"figurewindow_fontsize":				[8, float],
	"figurewindow_xlabel":					["", str],
	"figurewindow_ylabel":					["", str],
	"figurewindow_xlabelpadding":			[0, float],
	"figurewindow_ylabelpadding":			[0, float],
	"figurewindow_figkwargs":				[{}, dict],
	"figurewindow_expkwargs":				[{"color": "black", "linewidth": 0.7}, dict],
	"figurewindow_catkwargs":				[{"linewidth": 0.7}, dict],
	"figurewindow_linkwargs":				[{"marker": "*", "linewidth": 0.7, "s": 6}, dict],
	"figurewindow_border":					["Automatic", str],
	"figurewindow_leftborder":				[0.1, float],
	"figurewindow_rightborder":				[0.9, float],
	"figurewindow_topborder":				[0.9, float],
	"figurewindow_bottomborder":			[0.1, float],
	"figurewindow_annotation":				["Like main plot", str],
	"figurewindow_customannotation":		["", str],
	
	"calibratewindow_color":				["#ff0000", Color],
	"calibratewindow_derivative":			[0, int],
	"calibratewindow_lineshape":			["Gauss", str],
	"calibratewindow_querycat":				["", str],
	"calibratewindow_queryexp":				["", str],
	"calibratewindow_blendwidth":			[0.1, float],
	"calibratewindow_fitwidth":				[2, float],

	"files_exp":							[{}, dict],
	"files_cat":							[{}, dict],
	"files_lin":							[{}, dict],
}

# Generate with following command, quotes is str of quotes.txt
# json.dumps(quotes.split("\n\n"))
quotes_str = '["Erstens kommt es anders und zweitens als man denkt\\n- Ekrem Bora", "In der Theorie willst du was ver\\u00e4ndern, aber in der Praxis geht alles schief.\\nIn der Theorie wollen wir zu viel, geben aber viel zu wenig, das ist kein fairer Deal.\\n- Anis M. Y. Ferchichi & Jochen Burchard", "Und bist du unten, dr\\u00fccken sie dich noch ein St\\u00fcck tiefer\\nNoch ein St\\u00fcck tiefer\\nNoch ein St\\u00fcck tiefer\\nNoch ein St\\u00fcck tiefer\\n- Anis M. Y. Ferchichi", "Da ting goes skrrrrahh, pap, pap, ka-ka-ka\\nSkidiki-pap-pap, and a pu-pu-pudrrrr-boom\\nSkya, du-du-ku-ku-dun-dun\\nPoom, poom, you dun know\\n- Michael Dapaah", "It wasn\'t me\\n- Orville Richard Burrell", "We are going so fast\\nBut time so slow\\n- Theory of Relativity", "D\\u00f6p, d\\u00f6p, d\\u00f6p, d\\u00f6p, d\\u00f6d\\u00f6d\\u00f6d\\u00f6d\\u00f6p\\nD\\u00f6, d\\u00f6d\\u00f6d\\u00f6p, d\\u00f6p, d\\u00f6d\\u00f6d\\u00f6d\\u00f6d\\u00f6p, d\\u00f6p\\n- Hans Peter Geerdes", "Skibadee, skibadanger\\nI am the rearranger\\n- Hans Peter Geerdes", "Respect to the Man in the Ice Cream Van\\n- Hans Peter Geerdes", "Hyper Hyper\\n- Hans Peter Geerdes", "If we could only slow the time\\nWe would have forever every night\\n- Don Pepijn Schipper", "Meine Stra\\u00dfenpoesie l\\u00f6st die Chaostheorie\\n- Mousa Amouei", "Die Parabel sie steigt, und zwar exponentiell\\n- Mohamed El Moussaoui", "Chuba chuba chuba chuba chuba chuba chubby.\\nI don\'t have any lines to go right here, so chuby Teletubby\\n- Marshall Bruce Mathers III", "Two things are infinite: The universe and human stupidity;\\nand I\\u2018m not sure about the universe\\n- Albert E", "Physics is like sex: sure, it may give some practical results, but that\'s not why we do it.\\n- Richard P. Feynman", "I do not think you can name many great inventions that have been made by married men.\\n- Nikola Tesla", "Those who are not shocked when they first come across quantum theory cannot possibly have understood it.\\n- Niels Bohr", "We\\u2019re not free in what we do, because we\\u2019re not free in what we want.\\n- Jonas Kahnwald", "What we know is a drop. What we don\\u2019t know is an ocean.\\n- Isaac Newton", "Das ist der Sound f\\u00fcr die echten M\\u00e4nner, die das hier h\\u00f6ren, wenn sie Pressluft h\\u00e4mmern\\n- Tarek Ebene, Nico Seyfrid & Maxim Dr\\u00fcner", "I accept that\\n- Chuck Marstein", "Many of life\'s failures are people who did not realize how close they were to success when they gave up\\n- Thomas A. Edison", "I find that the harder I work, the more luck I seem to have\\n- Thomas Jefferson", "Whether you think you can or you think you can\'t, you\'re right\\n- Henry Ford", "Life is never fair, and perhaps it is a good thing for most of us that it is not\\n- Oscar Wilde", "Only a life lived for others is a life worthwhile\\n- Albert Einstein", "Imagination is more important than knowledge\\n- Albert Einstein", "My mama always said, \\u2018Life was like a box of chocolates. You never know what you\\u2019re gonna get\\n- Forrest", "Before you marry a person, you should first make them use a computer with slow Internet to see who they really are\\n- Will Ferrell", "I don\\u2019t believe in astrology; I\\u2019m a Sagittarius and we\\u2019re skeptical\\n- Arthur C. Clarke", "If you think you are too small to make a difference, try sleeping with a mosquito\\n- Dalai Lama", "People who think they know everything are a great annoyance to those of us who do\\n- Isaac Asimov", "They asked me how well I understood theoretical physics, I said I had a theoretical degree in physics.\\nThey said welcome aboard.\\n- Fantastic"]'

def start():
	global main
	main = Main()
	main.gui()

if __name__ == '__main__':
	start()
