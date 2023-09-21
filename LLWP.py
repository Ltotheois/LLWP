#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Luis Bonah
# Description: Loomis-Wood Plotting Software

CREDITSSTRING = """Made by Luis Bonah

As this programs GUI is based on PyQt5, which is GNU GPL v3 licensed, this program is also licensed under GNU GPL v3 (See the bottom paragraph).

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
PYQT_MAX = 2147483647

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

from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *

import matplotlib
from matplotlib import style, figure
from matplotlib.backends.backend_qt5agg import FigureCanvas, NavigationToolbar2QT

import warnings
warnings.simplefilter('ignore', np.RankWarning)

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
		rite = True
		q = main.plotwidget.working
		ind = main.plotwidget.indicator
		if q.empty():
			ind.setText("<span style='font-weight: 600;'>Working...</span>")
		q.put(1)
		try:
			res = func(self, *args, **kwargs)
		except Exception as E:
			res = str(E)
			rite = False
			raise
		finally:
			q.get()
			q.task_done()
			if q.empty():
				ind.setText("Ready")
			if rite:
				return(res)
	return(wrapper)

locks = {key: threading.RLock() for key in (
  "exp_df", "cat_df", "lin_df", "ser_df", "windows", "axs", "pipe", "currThread", "fitting", "peaks"
)}




##
## Main Class, Window and Widget
##
class Main():
	def __init__(self):
		self.open_windows = {}
		self.signalclass = SignalClass()
		self.config = Config(self.signalclass.updateconfig.emit)
		
		self.setup_dfs()
		self.loadoptions()
	
	def gui(self):
		sys.excepthook = except_hook
		self.app = QApplication(sys.argv)
		self.app.setStyle("Fusion")
		
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
			self.loadproject(sys.argv[1])
		self.mainwindow.show()
		
		if self.messages:
			main.notification("\n".join(self.messages))
		self.change_style(self.config["layout_theme"])
		sys.exit(self.app.exec_())
		
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
		if self.config["flag_alwaysshowlog"]:
			main.mainwindow.logwindow.setVisible(True)
			main.mainwindow.logwindow.raise_()
		self.signalclass.writelog.emit(output)

	def update_plot(self, dict_={}):
		self.signalclass.updateplot.emit(dict_)

	def loadoptions(self, fname=None):
		if not fname:
			self.config.update({key: value[0] for key, value in config_specs.items()})
			fname = f"{APP_TAG}.ini"
			
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
						value = class_(value)
						self.config[fullkey] = value
					except Exception as E:
						message = f"The value for the option {fullkey} from the option file was not understood."
						self.messages.append(message)
						print(message)
				else:
					self.config[fullkey] = value

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

			output_dict[category][name] = value

		del output_dict["Files"]

		config_parser = configparser.ConfigParser(interpolation=None)
		for section in output_dict:
			config_parser.add_section(section)
			for key in output_dict[section]:
				config_parser.set(section, key, str(output_dict[section][key]))

		with open(f"{APP_TAG}.ini", "w+", encoding="utf-8") as file:
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
				df.drop(df.index, inplace = True)
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
				self.load_file_core(fname, type, config_updates, results, errors, do_QNs)
		
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
			
			if len(fnames) != 0:
				xranges = [options["xrange"] for options in self.config["files_exp"].values()]
				xranges = sorted(xranges, key=lambda x: x[0])
				for i in range(len(xranges)-1):
					if xranges[i][1] > xranges[i+1][0]:
						self.notification(f"<span style='color:#eda711;'>WARNING</span>: Your spectra files are overlapping. This may cause zig-zag patterns. Think about using the Spectra Resolver module")
						break
		elif type == "cat":
			self.yrange_exp = np.array((df["y"].min(), df["y"].max()))

			if do_QNs and len(df) and self.config["flag_autosetqns"]:
				qn_labels = ['qnu1', 'qnu2', 'qnu3', 'qnu4', 'qnu5', 'qnu6']
				QNs = len(qn_labels)
				for i in range(len(qn_labels)):
					tmp_label = qn_labels[i]
					tmp_unique = df[tmp_label].unique()
					if len(tmp_unique) == 1 and tmp_unique[0] == np.iinfo(np.int64).min:
						QNs = i
						break
				
				self.config["series_qns"] = QNs
				self.notification(f"After analysing your cat files the number of QNs was set to {QNs}.")

		elif type == "lin":
			pass
				
		errors = list(errors.queue)
		self.signalclass.updatewindows.emit()
		if skip_update != True:
			self.plotwidget.set_data()
		if len(fnames) != 0:
			error_text = f"<span style='color:#ff0000;'>ERROR</span>: Reading {type.capitalize()} files not successful. " if len(errors) != 0 else ''
			self.notification(f"{error_text}Read {str(len(fnames)-len(errors))+'/' if len(errors) != 0 else ''}{len(fnames)} {type.capitalize()} files successfully.")
	
	def load_file_core(self, fname, type, config_updates, results, errors, do_QNs):
		try:
			if not os.path.isfile(fname):
				errors.put(fname)
				self.notification(f"<span style='color:#ff0000;'>ERROR</span>: The file {fname} could not be found. Please check the file.")
			if os.path.getsize(fname) == 0:
				return
			
			options = self.config[f"files_{type}"].get(fname, {})
			extension = os.path.splitext(fname)[1]
			
			if options.get("color") == None:
				options["color"] = self.config[f"color_{type}"]
			
			if type == "exp":
				args = (chr(self.config["flag_separator"]), self.config["flag_xcolumn"], self.config["flag_ycolumn"], options.get("invert", False), False)
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
							data[column] = np.iinfo(np.int64).min
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
							data[column] = np.iinfo(np.int64).min
					data = data[lin_dtypes.keys()]
				else:
					data = pyckett.lin_to_df(fname, False)

			config_updates.put({fname: options})
			results.put(data)
		except Exception as E:
			self.notification(f"<span style='color:#ff0000;'>ERROR</span>: There occurred an error when loading the {type.capitalize()} File {fname}. Please check the file.")
			if self.config["flag_debug"] == True:
				tb.print_exc()
			errors.put(fname)
			raise

	@threading_d
	def reread_files(self, do_QNs=False):
		kwargs = {"reread": True, "skip_update":True, "do_QNs": do_QNs}
		threads = []
		for type in ("exp", "cat", "lin"):
			threads.append(self.load_file(type, **kwargs))
		
		for thread in threads:
			thread.join()
		self.plotwidget.set_data()

	@synchronized_d(locks["lin_df"])
	def save_lines_lin(self, path = None, force_noappend=False, force_lin=False):
		append = self.config["flag_appendonsave"]
		options = {"options": QFileDialog.DontConfirmOverwrite} if append else {}

		if not path:
			path, ext = QFileDialog.getSaveFileName(None, 'Save file', '', **options)
			if not path:
				return
		
		catalogue = self.new_df
		output = []

		handle = "a+" if append and  not force_noappend else "w+"

		with open(path, handle, encoding="utf-8") as file:
			custom_format = self.config["flag_assignmentsavefmt"]
			if custom_format and not force_lin:
				np.savetxt(file, catalogue[custom_format["names"]], delimiter=custom_format.get("delimiter", " "), fmt=custom_format.get("format", '%.18e'))
			else:
				file.write(pyckett.df_to_lin(catalogue))
		self.notification(f"Newly assigned lines were saved to the file {path}.")

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
					for file in self.config[f"files_{type}"].keys():
						self.load_file(type, add_files=[file], keep_old=True)
				
		except Exception as E:
			self.notification(f"The command '{command}' failed with the Exception '{E}'.")
			raise

	def get_visible_data(self, type, xrange=None, binning=None, force_all=False):
		if type == "exp":
			with locks["exp_df"]:
				df = self.exp_df.copy()
				fd = self.config["files_exp"]
		elif type == "cat":
			with locks["cat_df"]:
				df = self.cat_df.copy()
				fd = self.config["files_cat"]
		elif type == "lin":
			with locks["lin_df"]:
				self.new_df["filename"] = "__lin_own_df__"
				df = pd.concat((self.lin_df, self.new_df), join="inner", ignore_index=True).sort_values("x")
				fd = self.config["files_lin"]

		if xrange != None:
			x_start = df["x"].searchsorted(xrange[0], side="left")
			x_stop  = df["x"].searchsorted(xrange[1], side="right")
			df = df.iloc[x_start:x_stop].copy()

		df_len = len(df.index)
		if df_len > main.config["plot_skipbinning"]:
			if binning != None and xrange != None:
				if binning == True:
					bins = self.config["plot_bins"]
				bin_width = (xrange[1]-xrange[0])/bins
				if bin_width == 0:
					bin_width = df_len

				df.loc[:,"binning"] = (df.loc[:,"x"]-xrange[0])//bin_width
				df = df.loc[df.sort_values(["y"]).drop_duplicates("binning", keep="last").sort_values(["x"]).index]

		if force_all != True:
			visible_files = {file for file in fd.keys() if not fd[file].get("hidden", False)}
			# Special Case Hide/Show catalogue files
			if type == "lin" and self.config["flag_hidecatalogue"] == False:
				visible_files.add("__lin_own_df__")

			if len(visible_files) != len(fd) + (type == "lin"):
				df.query("filename in @visible_files", inplace=True)
		return(df)

	def return_df(self, type):
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

	def change_style(self, style=None):
		styles = ["light", "dark", "custom"]
		app = self.app
		if style == None:
			self.config["layout_theme"] = styles[(styles.index(self.config["layout_theme"])+1)%len(styles)]
		elif style in styles:
			self.config["layout_theme"] = style
		else:
			self.config["layout_theme"] = styles[0]
		
		if self.config["layout_owntheme"] == {} and self.config["layout_theme"] == "custom":
			self.config["layout_theme"] = "light"

		if self.config["layout_theme"] == "light":
			palette = app.style().standardPalette()
			mplstyles = ("default", "white")

		elif self.config["layout_theme"] == "dark" or self.config["layout_theme"] == "custom":
			colors = {
				"window":				QColor(53, 53, 53),
				"windowText":			QColor(255, 255, 255),
				"base":					QColor(35, 35, 35),
				"alternateBase":		QColor(53, 53, 53),
				"toolTipBase":			QColor(25, 25, 25),
				"toolTipText":			QColor(255, 255, 255),
				"placeholderText":		QColor(100, 100, 100),
				"text":					QColor(255, 255, 255),
				"button":				QColor(53, 53, 53),
				"buttonText":			QColor(255, 255, 255),
				"brightText":			Qt.red,
				"light":				QColor(255, 255, 255),
				"midlight":				QColor(200, 200, 200),
				"mid":					QColor(150, 150, 150),
				"dark":					QColor(50, 50, 50),
				"shadow":				QColor(0, 0, 0),
				"highlight":			QColor(42, 130, 218),
				"highlightedText":		 QColor(35, 35, 35),
				"link":					QColor(42, 130, 218),
				"linkVisited":			QColor(42, 130, 218),
				
				"disabledButtonText":	Qt.darkGray,
				"disabledWindowText":	Qt.darkGray,
				"disabledText":			Qt.darkGray,
				"disabledLight":		QColor(53, 53, 53),
				
				"mplstyles":			("dark_background", "black"),
			}
			
			if self.config["layout_theme"] == "custom":
				colors.update(self.config["layout_owntheme"])
			
			tmp_dict = {
				"window":				(QPalette.Window,),
				"windowText":			(QPalette.WindowText,),
				"base":					(QPalette.Base,),
				"alternateBase":		(QPalette.AlternateBase,),
				"toolTipBase":			(QPalette.ToolTipBase,),
				"toolTipText":			(QPalette.ToolTipText,),
				"placeholderText":		(QPalette.PlaceholderText,),
				"text":					(QPalette.Text,),
				"button":				(QPalette.Button,),
				"buttonText":			(QPalette.ButtonText,),
				"brightText":			(QPalette.BrightText,),
				"light":				(QPalette.Light,),
				"midlight":				(QPalette.Midlight,),
				"dark":					(QPalette.Dark,),
				"mid":					(QPalette.Mid,),
				"shadow":				(QPalette.Shadow,),
				"highlight":			(QPalette.Highlight,),
				"highlightedText":		(QPalette.HighlightedText,),
				"link":					(QPalette.Link,),
				"linkVisited":			(QPalette.LinkVisited,),
				
				"disabledButtonText":	(QPalette.Disabled, QPalette.ButtonText),
				"disabledWindowText":	(QPalette.Disabled, QPalette.WindowText),
				"disabledText":			(QPalette.Disabled, QPalette.Text),
				"disabledLight":		(QPalette.Disabled, QPalette.Light),
			}
			
			mplstyles = colors["mplstyles"]
			palette = QPalette()
			for key, values in tmp_dict.items():
				palette.setColor(*values, colors[key])
			
		
		app.setPalette(palette)
		matplotlib.style.use(mplstyles[0])
		matplotlib.rc('font', **self.config["plot_font_dict"])
		main.config.register("plot_font_dict", lambda: matplotlib.rc('font', **self.config["plot_font_dict"]))
		self.plotwidget.fig.patch.set_facecolor(mplstyles[1])
		self.plotwidget.create_plots()

	@synchronized_d(locks["windows"])
	def open_window(self, window, *args):
		if window not in self.open_windows:
			self.open_windows[window] = window(window.__name__.lower(), *args)
		self.open_windows[window].show()
		self.open_windows[window].activateWindow()
		return(self.open_windows[window])

	@synchronized_d(locks["lin_df"])
	def catalogue_table_delete(self, all=False):
		table = main.mainwindow.cataloguewindow.catalogueTable
		model = main.mainwindow.cataloguewindow.catalogueModel
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
		self.setFocusPolicy(Qt.StrongFocus)
		self.setWindowTitle(APP_TAG)
		self.setAcceptDrops(True)
		self.shortcuts()

		geometry = main.config.get("windowgeometry_mainwindow")
		if geometry:
			self.setGeometry(*json.loads(geometry))
		
		try:
			main.app.setWindowIcon(QIcon(f"{APP_TAG}.svg"))
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
		self.cataloguewindow = CatalogueWindow()
		self.logwindow = LogWindow()
		self.quotewindow = QuoteWindow()

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
		toggleaction_referenceseries = self.referenceserieswindow.toggleViewAction()
		toggleaction_referenceseries.setShortcut("Shift+3")
		toggleaction_catalogue = self.cataloguewindow.toggleViewAction()
		toggleaction_catalogue.setShortcut("Shift+4")
		toggleaction_log = self.logwindow.toggleViewAction()
		toggleaction_log.setShortcut("Shift+5")

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
				QQ(QAction, parent=self, text="&Edit Spectrum Files", shortcut="Shift+6", tooltip="See overview of current spectrum files and their options", change=lambda x: main.open_window(ExpWindow)),
				QQ(QAction, parent=self, text="&Edit Cat Files", shortcut="Shift+7", tooltip="See overview of current .cat files and their options", change=lambda x: main.open_window(CatWindow)),
				QQ(QAction, parent=self, text="&Edit Lin Files", shortcut="Shift+8", tooltip="See overview of current .lin files and their options", change=lambda x: main.open_window(LinWindow)),
				None,
				QQ(QAction, parent=self, text="&Save current values as default", shortcut="Ctrl+D", tooltip="Save current configuration as default", change=lambda x: main.saveoptions()),
				None,
				QQ(QAction, parent=self, text="&Save Project", change=lambda x: main.saveproject(), shortcut="Ctrl+S"),
				QQ(QAction, parent=self, text="&Load Project", change=lambda x: main.loadproject(), shortcut="Ctrl+O"),
				None,
				QQ(QAction, parent=self, text="&Quit", change=self.close),
			),
			"View": (
				QQ(QAction, parent=self, text="&Change Style", tooltip="Change between light, dark and custom theme", change=lambda x: main.change_style()),
				None,
				QQ(QAction, "layout_mpltoolbar", parent=self, text="&Show MPL Toolbar", shortcut="Shift+1", tooltip="Show or hide toolbar to edit or save the plot canvas", checkable=True),
				toggleaction_configureplots,
				toggleaction_referenceseries,
				toggleaction_catalogue,
				toggleaction_log,
				None,
				QQ(QAction, "layout_limitannotationtosingleline", parent=self, text="&Annotate is Single Line", tooltip="Switch between the onhover annotation line being one line (for consistent width and height) or multiple lines (if many lines are at the same position and you want to see all of them)", checkable=True),
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
				QQ(QAction, parent=self, text="&Manual Assign", tooltip="Manually assign and add line to the assigned lines table, useful when always assign is not selected", change=lambda x: self.plotwidget.manual_assign(), shortcut="Ctrl+Return"),
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
				QQ(QAction, parent=self, text="&Show Seriesfinder", shortcut="Ctrl+3", tooltip="Press to open the seriesfinder module, which allows to explore the .cat lines, e.g. finding the most intense transitions", change=lambda x: main.open_window(SeriesfinderWindow)),
				QQ(QAction, parent=self, text="&Peakfinder", shortcut="Ctrl+4", tooltip="Press to open the peakfinder module, which allows to find the peaks of your experimental spectrum and check, which peaks are not assigned yet", change=lambda x: main.open_window(PeakfinderWindow)),
				QQ(QAction, parent=self, text="&Residuals", shortcut="Ctrl+5", tooltip="Press to open the Residuals module", change=lambda x: main.open_window(ResidualsWindow)),
				QQ(QAction, parent=self, text="&Energy Levels", shortcut="Ctrl+6", tooltip="Press to open the Energy Levels module", change=lambda x: main.open_window(EnergyLevelsWindow)),
				QQ(QAction, parent=self, text="&Spectra Resolver", shortcut="Ctrl+7", tooltip="Press to open the Spectra Resolver module", change=lambda x: main.open_window(SpectraResolverWindow)),
				QQ(QAction, parent=self, text="&Show Lineshape", shortcut="Ctrl+8", tooltip="Press to open the lineshape module, which allows to simulate different line profiles from the .cat stick spectrum", change=lambda x: main.open_window(LineshapeWindow)),
				QQ(QAction, parent=self, text="&Series Fit", shortcut="Ctrl+9", tooltip="Open Series Fit Window, which allows to try out different combinations of transitions as a series", change=lambda x: main.open_window(SeriesFitWindow)),
				QQ(QAction, parent=self, text="&Create Report", tooltip="Open Report Window, which allows to summarise your analysis", change=lambda x: main.open_window(ReportWindow)),
				QQ(QAction, parent=self, text="&Create Figure", tooltip="Create publication quality figure", change=lambda x: main.open_window(FigureWindow)),
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

			extension = os.path.splitext(file)[1]
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
			t = main.loadproject(project)
			t.join()
			
		main.plotwidget.set_data()

	def closeEvent(self, *args, **kwargs):
		if not main.new_df.empty:
			main.save_lines_lin(f"{APP_TAG}.lin", force_noappend=True, force_lin=True)
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
		main.config.register("plot_dpi", lambda: self.fig.set_dpi(main.config["plot_dpi"]))
		self.plotcanvas = FigureCanvas(self.fig)
		self.plotcanvas.setMinimumHeight(200)
		self.plotcanvas.setMinimumWidth(200)
		self.fig.canvas.mpl_connect("motion_notify_event", self.on_hover)
		self.cache_positions = None

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

		self.toplabel = QQ(QLabel, text="", wordwrap=True, maxHeight=(12 if main.config["layout_limitannotationtosingleline"] else 10000))
		self.indicator = QQ(QLabel, text="Ready", textFormat=Qt.RichText)
		self.working = queue.Queue()

		main.config.register("layout_limitannotationtosingleline", lambda: self.toplabel.setMaximumHeight(12 if main.config["layout_limitannotationtosingleline"] else 16777215))

		toplayout.addWidget(self.toplabel, 1)
		toplayout.addWidget(self.indicator)
		
		rows_cols_elements = (
			QQ(QLabel, text="||  Plots: ", visible=main.config["flag_showmainplotrowscols"]),
			QQ(QSpinBox, "plot_rows", range=(1, PYQT_MAX), maxWidth=45, visible=main.config["flag_showmainplotrowscols"]),
			QQ(QLabel, text="x", visible=main.config["flag_showmainplotrowscols"]),
			QQ(QSpinBox, "plot_cols", range=(1, PYQT_MAX), maxWidth=45, visible=main.config["flag_showmainplotrowscols"]),
		)
		for elem in rows_cols_elements:
			toplayout.addWidget(elem)
			main.config.register("flag_showmainplotrowscols", lambda elem=elem: elem.setVisible(main.config["flag_showmainplotrowscols"]))

		width_elements = (
			QQ(QLabel, text="    Width: ", visible=main.config["flag_showmainplotwidth"]),
			QQ(QDoubleSpinBox, "plot_width", range=(0, PYQT_MAX), maxWidth=75, visible=main.config["flag_showmainplotwidth"], change=lambda: main.config["plot_coupled"] and self.set_data()),
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
		main.config.register(("series_annotate_xs", "plot_annotate", "plot_bins"), lambda: self.set_data())
	
	def gui(self):
		self.create_plots().join()
		main.config.register(("plot_rows", "plot_cols"), self.create_plots)
		main.config.register("series_references", self.set_data)

	def set_offset(self, value):
		if main.config["plot_coupled"]:
			offset = main.config["plot_offset"]
			width = main.config["plot_width"]
		else:
			offset = self.axs["offset"][self.lcp]
			width = self.axs["width"][self.lcp]
		
		if value == "+":
			offset += width/4
		elif value == "-":
			offset -= width/4
		elif value == "++":
			offset += width/2
		elif value == "--":
			offset -= width/2
		else:
			offset = value
		
		if main.config["plot_coupled"]:
			main.config["plot_offset"] = offset
		else:
			self.axs["offset"][self.lcp] = offset
		self.set_data()

	def get_offset(self, index=None):
		if main.config["plot_coupled"]:
			return main.config["plot_offset"]
		else:
			if index is None:
				index = self.get_current_plot()
			return self.axs["offset"][index]
	
	def reset_offsets(self):
		shape = self.axs["ax"].shape
		self.axs["offset"] = np.zeros(shape)
		self.axs["width"]  = np.full(shape, main.config["plot_width"], dtype=np.float64)
		main.config["plot_offset"] = 0
		self.set_data()
	
	def set_width(self, value, absolute=True):
		if main.config["plot_coupled"]:
			width = main.config["plot_width"]
		else:
			width = self.axs["width"][self.lcp]
		
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
		
		if main.config["plot_coupled"]:
			main.config["plot_width"] = width
		else:
			self.axs["width"][self.lcp] = width
		self.set_data()
	
	def get_widths(self):
		if main.config["plot_coupled"]:
			return np.full(self.axs["ax"].shape, main.config["plot_width"])
		else:
			return self.axs["width"]

	def get_positions(self, return_qns=False):
		shape = self.axs["ax"].shape
		positions = np.zeros(shape)
		if return_qns:
			allqns = np.empty(shape, dtype=object)
		references = main.config["series_references"]
		
		for i, reference in zip(range(shape[1]), references):
			method = reference["method"]
			if method == "Transition":
				qns = self.get_qns(reference["transition"])
				if return_qns:
					allqns[:, i] = qns
				file = reference["transition"].get("file")
				positions[:, i] = self.get_position_from_qns(qns, file)
			elif method == "List":
				if return_qns:
					# qns = reference["list"]["qns"]
					allqns[:, i] = None
				i0 = reference["list"]["i0"]
				xs_all = reference["list"]["xs"]
				
				xs = np.zeros(shape[0])
				if xs_all:
					imax = min(len(xs), len(xs_all)-i0)
					xs[:imax] = xs_all[i0:imax+i0]
				positions[:, i] = xs[::-1]
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
			
			self.cache_positions = positions
		if return_qns:
			return(positions, allqns)
		else:
			return(positions)
	
	def get_seriesmethods(self):
		shape = self.axs["ax"].shape
		references = main.config["series_references"]
		seriesmethods = [reference["method"] for reference in references]
		return(seriesmethods)
	
	def get_qns(self, transition):
		nop = self.axs["ax"].shape[0]
		noq = main.config["series_qns"]
		
		qnus, qnls = np.array(transition["qnus"][:noq]), np.array(transition["qnls"][:noq])
		diffs = np.array(transition["qndiffs"][:noq]) if transition["increase_freely"] else np.array(transition["qnincrs"][:noq])
		
		qns = []
		
		for i in range(nop):
			ind = nop-i-1
			upper, lower = qnus+ind*diffs, qnls+ind*diffs
			qns.append((upper, lower))
		return(qns)
	
	def get_position_from_qns(self, qns, file=None):
		positions = []
		for qnus, qnls in qns:
			cond_upper = [f"(qnu{i+1} == {qn})" for i, qn in enumerate(qnus)]
			cond_lower = [f"(qnl{i+1} == {qn})" for i, qn in enumerate(qnls)]
			condition  = " & ".join(cond_upper + cond_lower)
			if file:
				condition += f" & (filename == @file)"
			
			vals = main.cat_df.query(condition)["x"].to_numpy()
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
						self.axs["annotation"][i, j] = None
			return
		
		for i in range(self.axs["ax"].shape[0]):
			for j in range(self.axs["ax"].shape[1]):
				annotation = self.axs["annotation"][i, j]
				x = xpos[i, j]
				qns = allqns[i, j]
				color = matplotlib.rcParams['text.color']
				
				if main.config["series_annotate_xs"] or qns is None:
					text = f"{{:{main.config['series_annotate_fmt']}}}".format(x)
				else:
					text = f"{', '.join([str(qn) if qn != np.iinfo(np.int64).min else '-' for qn in qns[0]])} ← {', '.join([str(qn) if qn != np.iinfo(np.int64).min else '-' for qn in qns[1]])}"
					
					lin_df = main.get_visible_data("lin")
					if len(lin_df.query(" and ".join([f"qnu{i+1} == {qn}" for i, qn in enumerate(qns[0])] + [f"qnl{i+1} == {qn}" for i, qn in enumerate(qns[1])]))):
						color = main.config["color_lin"]
			
				if not annotation:
					ax  = self.axs["ax"][i, j]
					self.axs["annotation"][i, j] = ax.text(**ann_dict, s=text, transform=ax.transAxes)
				else:
					self.axs["annotation"][i, j].set_text(text)
				self.axs["annotation"][i, j].set_color(color)

	def check_blends(self, index, dict_):
		if main.config["series_blenddialog"] and self.get_seriesmethods()[index[1]]:
			blendwidth = main.config["series_blendwidth"]
			xrange = (dict_["xpre"]-blendwidth/2, dict_["xpre"]+blendwidth/2)
			entries = main.get_visible_data("cat", xrange=xrange)
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
		
		if all([main.config["plot_hover"], x, y, event.inaxes]):
			text = f"({x=:.2f}, {y=:.2f}) "

			cutoff = main.config["plot_hover_cutoff"]
			df = main.get_visible_data("cat", xrange=(x-cutoff, x+cutoff))
			if len(df):
				df["dist"] = abs(df["x"] - x)
				df.sort_values("dist", inplace=True)
				
				transitions = []
				noq = main.config["series_qns"]
				maxdist = df.iloc[0]["dist"]
				for i, row in df.iterrows():
					if maxdist < row["dist"]:
						break
					qnus = [row[f"qnu{i+1}"] for i in range(noq)]
					qnls = [row[f"qnl{i+1}"] for i in range(noq)]
					transitions.append(f"{', '.join(map(str, map(int, qnus)))} ← {', '.join(map(str, map(int, qnls)))}")
					
				text += " || ".join(transitions)
		else:
			text = ""
		self.toplabel.setText(text)

	@synchronized_d(locks["axs"])
	def on_range(self, xmin, xmax, index):
		self.lcp = index
		axrange = self.axs["ax"][index].get_xlim()
		if xmax == xmin or xmax > axrange[1] or xmin < axrange[0]:
			return
		
		shift = (QApplication.keyboardModifiers() == Qt.ShiftModifier)
		if shift or main.config["fit_alwaysfit"]:
			xmiddle, xuncert = self.fit_peak(xmin, xmax, index)
			xpre = self.cache_positions[index]
			
			dict_ = {
				"x":		xmiddle,
				"error":	xuncert,
				"xpre":		xpre
			}
			reference = main.config["series_references"][index[1]]
			if reference["method"] == "Transition":
				qns = self.get_qns(reference["transition"])[index[0]]
				for i, qnu, qnl in zip(range(6), qns[0], qns[1]):
					dict_[f"qnu{i+1}"] = qnu
					dict_[f"qnl{i+1}"] = qnl
	
			if not main.config["fit_alwaysassign"]:
				self.lastassignment = dict_
			elif self.check_blends(index, dict_):
				return
			else:
				self.assign(index, dict_)
		else:
			self.set_offset( (xmin+xmax)/2 - self.cache_positions[index])
			self.set_width(xmax-xmin)
	
	def fit_peak(self, xmin, xmax, index=None):
		data = main.get_visible_data("exp", xrange=(xmin, xmax))
		exp_xs = data["x"].to_numpy()
		exp_ys = data["y"].to_numpy()
		fit_xs = np.linspace(xmin, xmax, main.config["fit_xpoints"])
		
		if len(exp_xs) < 2:
			main.notification("<span style='color:#eda711;'>WARNING</span>: You did select less than two points of your spectrum, this fit will not work.")
			raise CustomError("Too few points selected")
		
		if self.fitline != None:
			self.fitline.remove()
			self.fitline = None
		if self.fitcurve != None:
			self.fitcurve.remove()
			self.fitcurve = None
		
		try:
			fitfunction = main.config["fit_function"]
			if fitfunction == "Pgopher":
				ymin, ymax = np.min(exp_ys), np.max(exp_ys)
				cutoff = ymax - (ymax-ymin)/2
				
				mask = (exp_ys >= cutoff)
				fit_xs = exp_xs[mask]
				fit_ys = exp_ys[mask]
				
				xmiddle = np.sum(fit_xs*fit_ys)/np.sum(fit_ys)
				xuncert = 0
			elif fitfunction == "Polynom":
				try:
					popt = np.polyfit(exp_xs, exp_ys, main.config["fit_polynomrank"])
				except Exception as E:
					popt = np.polyfit(exp_xs, exp_ys, main.config["fit_polynomrank"])
				polynom = np.poly1d(popt)
				fit_ys = polynom(fit_xs)
				
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

				xmiddle = fit_xs[np.argmax(fit_ys)]
				xuncert = 0 # @Luis: Find real error for this part
			
			else:
				x0 = (xmin+xmax)/2
				ymin, ymax = np.min(exp_ys), np.max(exp_ys)
				ymean = np.mean(exp_ys)
				y0 = ymax-ymin
				w0 = (xmax-xmin)
				
				if main.config["fit_offset"]:
					function, p0, bounds = {
						"Gauss":				(lambda *x: lineshape("Gauss", 0, *x[:-1])+x[-1],   (x0, y0, w0, ymean),      ((xmin, 0, 0, ymin),    (xmax, 3*y0, 10*w0, ymax))),
						"Lorentz":				(lambda *x: lineshape("Lorentz", 0, *x[:-1])+x[-1], (x0, y0, w0, ymean),      ((xmin, 0, 0, ymin),    (xmax, 3*y0, 10*w0, ymax))),
						"Voigt":				(lambda *x: lineshape("Voigt", 0, *x[:-1])+x[-1],   (x0, y0, w0, w0, ymean),  ((xmin, 0, 0, 0, ymin), (xmax, 3*y0, 5*w0, 5*w0, ymax))),
						"Gauss 1st Derivative":		(lambda *x: lineshape("Gauss", 1, *x[:-1])+x[-1],   (x0, y0, w0, ymean),      ((xmin, 0, 0, ymin),    (xmax, 3*y0, 10*w0, ymax))),
						"Lorentz 1st Derivative":		(lambda *x: lineshape("Lorentz", 1, *x[:-1])+x[-1], (x0, y0, w0, ymean),      ((xmin, 0, 0, ymin),    (xmax, 3*y0, 10*w0, ymax))),
						"Voigt 1st Derivative":		(lambda *x: lineshape("Voigt", 1, *x[:-1])+x[-1],   (x0, y0, w0, w0, ymean),  ((xmin, 0, 0, 0, ymin), (xmax, 3*y0, 5*w0, 5*w0, ymax))),
						"Gauss 2nd Derivative":		(lambda *x: lineshape("Gauss", 2, *x[:-1])+x[-1],   (x0, y0, w0, ymean),      ((xmin, 0, 0, ymin),    (xmax, 3*y0, 10*w0, ymax))),
						"Lorentz 2nd Derivative":		(lambda *x: lineshape("Lorentz", 2, *x[:-1])+x[-1], (x0, y0, w0, ymean),      ((xmin, 0, 0, ymin),    (xmax, 3*y0, 10*w0, ymax))),
						"Voigt 2nd Derivative":		(lambda *x: lineshape("Voigt", 2, *x[:-1])+x[-1],   (x0, y0, w0, w0, ymean),  ((xmin, 0, 0, 0, ymin), (xmax, 3*y0, 5*w0, 5*w0, ymax))),
					}.get(fitfunction)
				else:
					function, p0, bounds = {
						"Gauss":				(lambda *x: lineshape("Gauss", 0, *x),   (x0, y0, w0),      ((xmin, 0, 0),    (xmax, 3*y0, 10*w0))),
						"Lorentz":				(lambda *x: lineshape("Lorentz", 0, *x), (x0, y0, w0),      ((xmin, 0, 0),    (xmax, 3*y0, 10*w0))),
						"Voigt":				(lambda *x: lineshape("Voigt", 0, *x),   (x0, y0, w0, w0),  ((xmin, 0, 0, 0), (xmax, 3*y0, 5*w0, 5*w0))),
						"Gauss 1st Derivative":		(lambda *x: lineshape("Gauss", 1, *x),   (x0, y0, w0),      ((xmin, 0, 0),    (xmax, 3*y0, 10*w0))),
						"Lorentz 1st Derivative":		(lambda *x: lineshape("Lorentz", 1, *x), (x0, y0, w0),      ((xmin, 0, 0),    (xmax, 3*y0, 10*w0))),
						"Voigt 1st Derivative":		(lambda *x: lineshape("Voigt", 1, *x),   (x0, y0, w0, w0),  ((xmin, 0, 0, 0), (xmax, 3*y0, 5*w0, 5*w0))),
						"Gauss 2nd Derivative":		(lambda *x: lineshape("Gauss", 2, *x),   (x0, y0, w0),      ((xmin, 0, 0),    (xmax, 3*y0, 10*w0))),
						"Lorentz 2nd Derivative":		(lambda *x: lineshape("Lorentz", 2, *x), (x0, y0, w0),      ((xmin, 0, 0),    (xmax, 3*y0, 10*w0))),
						"Voigt 2nd Derivative":		(lambda *x: lineshape("Voigt", 2, *x),   (x0, y0, w0, w0),  ((xmin, 0, 0, 0), (xmax, 3*y0, 5*w0, 5*w0))),
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
			"qnu1":			np.iinfo(np.int64).min,
			"qnu2":			np.iinfo(np.int64).min,
			"qnu3":			np.iinfo(np.int64).min,
			"qnu4":			np.iinfo(np.int64).min,
			"qnu5":			np.iinfo(np.int64).min,
			"qnu6":			np.iinfo(np.int64).min,
			"qnl1":			np.iinfo(np.int64).min,
			"qnl2":			np.iinfo(np.int64).min,
			"qnl3":			np.iinfo(np.int64).min,
			"qnl4":			np.iinfo(np.int64).min,
			"qnl5":			np.iinfo(np.int64).min,
			"qnl6":			np.iinfo(np.int64).min,
			"x":			np.iinfo(np.int64).min,
			"error":		0,
			"weight":		1,
			"comment":		"",
			"filename":		"__lin_own_df__",
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
			if rc == True:
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
			if annotation:
				annotation.set_color(main.config["color_lin"])
				ax.draw_artist(annotation)
			
			lin_plot = self.axs["lin_plot"][index]
			if lin_plot:
				offsets = lin_plot.get_offsets()
				offsets = np.concatenate([offsets, np.array((lin_dict["x"], 0) ,ndmin=2)])
				lin_plot.set_offsets(offsets)
				
				ax.draw_artist(lin_plot)
			with locks["axs"]:
				self.plotcanvas.update()
				self.plotcanvas.flush_events()

	def get_current_plot(self):
		lcp = self.lcp
		shape = self.axs["ax"].shape
		if 0 <= lcp[0] < shape[0] and 0 <= lcp[1] < shape[1]:
			return(lcp)
		else:
			lcp = (0, 0)
			return(lcp)
	
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

			breakpoint(ownid, self.create_plots_id)
			
			tmp = self.fig.subplots(rows, cols, gridspec_kw=main.config["plot_matplotlibkwargs"], squeeze=False)
			
			breakpoint(ownid, self.create_plots_id)
			
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
			
			breakpoint(ownid, self.create_plots_id)
			
			for i, row in enumerate(tmp):
				for j, ax in enumerate(row):
					self.axs["ax"][i, j] = ax
					self.axs["lin_plot"][i, j] = ax.scatter([], [], color=main.config["color_lin"], marker="*")
					self.axs["span"][i, j] = matplotlib.widgets.SpanSelector(ax, lambda xmin, xmax, i=i, j=j:self.on_range(xmin, xmax, (i, j)), 'horizontal')
					ax.yaxis.set_visible(False)
					ax.xaxis.set_visible(False)
			
			breakpoint(ownid, self.create_plots_id)
			
			for bottomax in self.axs["ax"][-1, :]:
				bottomax.xaxis.set_visible(True)
			
			breakpoint(ownid, self.create_plots_id)
			
			main.signalclass.createdplots.emit()
			self.set_data()
		except CustomError as E:
			pass

	def manual_draw(self):
		self.set_data().join()
		if not main.config["flag_automatic_draw"]:
			with locks["axs"]:
				self.plotcanvas.draw()

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
			# set x-ranges
			# @Luis: Think about optimizations
			xpos, qns = self.get_positions(return_qns=True)
			widths = self.get_widths()
			
			xmin = np.zeros(xpos.shape)
			xmax = np.zeros(xpos.shape)
			
			for i in range(self.axs["ax"].shape[0]):
				for j in range(self.axs["ax"].shape[1]):
					ax = self.axs["ax"][i, j]
					offset = self.get_offset((i, j))
					x, width = xpos[i, j]+offset, widths[i, j]/2
					xrange = (x-width, x+width)
					ax.set_xlim(*xrange)
					
					xmin[i, j], xmax[i, j] = xrange
			
			# set ticks for bottom row
			i = -1
			for j in range(self.axs["ax"].shape[1]):
				ax = self.axs["ax"][i, j]
				ticks = np.linspace(xmin[i, j], xmax[i, j], main.config["plot_ticks"])
				if main.config["plot_offsetticks"] == 1:
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
				

			breakpoint(ownid, self.set_data_id)
			
			# set data and set y-ranges of data
			bins = main.config["plot_bins"]
			scaling = main.config["plot_yscale"]
			
			dataframes = [main.get_visible_data(x) for x in ("exp", "cat", "lin")]
			files_dicts = [main.config[x] for x in ("files_exp", "files_cat", "files_lin")]

			for i in range(self.axs["ax"].shape[0]):
				for j in range(self.axs["ax"].shape[1]):
					breakpoint(ownid, self.set_data_id)
					ax = self.axs["ax"][i, j]

					for datatype, dataframe, files in zip(("exp", "cat", "lin"), dataframes, files_dicts):
						minindex = dataframe["x"].searchsorted(xmin[i, j], side="left")
						maxindex = dataframe["x"].searchsorted(xmax[i, j], side="right")
						
						dataframe = dataframe.iloc[minindex:maxindex].copy()

						tmp_length = len(dataframe.index)
						if tmp_length > bins:
							tmp_binning = (xmax[i, j]-xmin[i, j])/bins

							if tmp_binning == 0:
								tmp_binning = tmp_length

							dataframe.loc[:,"binning"] = (dataframe.loc[:,"x"]-xmin[i, j])//tmp_binning
							dataframe = dataframe.loc[dataframe.sort_values(["y"]).drop_duplicates("binning", keep="last").sort_values(["x"]).index]

						xs = dataframe["x"].to_numpy()
						ys = dataframe["y"].to_numpy() if datatype != "lin" else 0*xs
						
						if datatype == "exp":
							if scaling == "Per Plot":
								if len(dataframe):
									yrange_exp = [dataframe["y"].min(), dataframe["y"].max()]
								else:
									yrange_exp = [-1, 1]
							segs = (((xs[i], ys[i]),(xs[i+1], ys[1+i])) for i in range(len(xs)-1))
							colors = create_colors(dataframe, files)
							coll = matplotlib.collections.LineCollection(segs, colors=colors)
							if self.axs["exp_plot"][i, j]:
								self.axs["exp_plot"][i, j].remove()
							self.axs["exp_plot"][i, j] = ax.add_collection(coll)
						elif datatype == "cat":
							if scaling == "Per Plot":
								if len(dataframe):
									yrange_cat = [dataframe["y"].min(), dataframe["y"].max()]
								else:
									yrange_cat = [-1, 1]
								ys = ys*yrange_exp[1]/yrange_cat[1]
							elif scaling in ["Global", "Custom"]:
								ys = ys*main.config["plot_expcat_factor"]*10**main.config["plot_expcat_exponent"]
							segs = (((xs[i], 0),(xs[i], ys[i])) for i in range(len(xs)))
							colors = create_colors(dataframe, files, xpos[i, j])
							coll = matplotlib.collections.LineCollection(segs, colors=colors)
							if self.axs["cat_plot"][i, j]:
								self.axs["cat_plot"][i, j].remove()
							self.axs["cat_plot"][i, j] = ax.add_collection(coll)
						elif datatype == "lin":
							tuples = list(zip(xs,ys))
							tuples = tuples if len(tuples)!=0 else [[None,None]]
							self.axs["lin_plot"][i, j].set_offsets(tuples)

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
			
			breakpoint(ownid, self.set_data_id)
			self.plot_annotations(xpos, qns)
			
			breakpoint(ownid, self.set_data_id)
			if main.config["flag_automatic_draw"]:
				with locks["axs"]:
					self.plotcanvas.draw()
		except CustomError as E:
			pass

	def wheelEvent(self, event):
		steps = event.angleDelta().y() // 120
		self.set_width(2 ** -steps, absolute=False)

	def change_fitcolor(self):
		color = QColorDialog.getColor(initial=QColor(rgbt_to_trgb(main.config["color_fit"])), options=QColorDialog.ShowAlphaChannel)
		color = Color(trgb_to_rgbt(color.name(QColor.HexArgb)))
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
		
		index = self.get_current_plot()
		xpos  = self.cache_positions[index]
		
		self.set_offset(xnew - xpos)
		
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
		
		main.mainwindow.addDockWidget(2, self)
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
		binningfield = QQ(QSpinBox, "plot_bins", range=(1, PYQT_MAX))
		[hbox1.addWidget(binninglabel), hbox1.addWidget(binningfield), hbox1.addStretch(1)]

		checkbox_coupled  = QQ(QCheckBox, "plot_coupled", text="Plots are coupled")
		checkbox_annotate = QQ(QCheckBox, "plot_annotate", text="Annotate plots")
		checkbox_hover    = QQ(QCheckBox, "plot_hover", text="Show nearest transition")

		
		[layout.addLayout(hbox1), layout.addItem(QSpacerItem(5, 5, QSizePolicy.Minimum, QSizePolicy.Expanding))]
		[layout.addWidget(widget) for widget in (checkbox_coupled, checkbox_annotate, checkbox_hover, )]

		checkbox_scale = QQ(QComboBox, "plot_yscale", options=("Per Plot", "Global", "Custom"), change=lambda x: main.signalclass.updateplot.emit())
		spinbox_scalemin = QQ(QDoubleSpinBox, "plot_yscale_min", range=(-PYQT_MAX, PYQT_MAX), maxWidth=80, change=lambda x: main.signalclass.updateplot.emit())
		spinbox_scalemax = QQ(QDoubleSpinBox, "plot_yscale_max", range=(-PYQT_MAX, PYQT_MAX), maxWidth=80, change=lambda x: main.signalclass.updateplot.emit())
		spinbox_scalecatfac = QQ(QDoubleSpinBox, "plot_expcat_factor", range=(0, PYQT_MAX), maxWidth=80, change=lambda x: main.signalclass.updateplot.emit())
		spinbox_scalecatexp = QQ(QSpinBox, "plot_expcat_exponent", range=(-PYQT_MAX, PYQT_MAX), prefix="*10^", maxWidth=80, change=lambda x: main.signalclass.updateplot.emit())

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
		grid.addWidget(spinbox_scalecatexp, 1, 2)
		grid.addWidget(scaleMinLabel, 2, 0)
		grid.addWidget(spinbox_scalemin, 2, 1)
		grid.addWidget(scaleMaxLabel, 2, 2)
		grid.addWidget(spinbox_scalemax, 2, 3)
		grid.setColumnStretch(7, 10)
		grid.setRowStretch(2, 10)
		layout.addItem(QSpacerItem(5, 5, QSizePolicy.Minimum, QSizePolicy.Expanding))
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
		
		self.tabs.setCornerWidget(QQ(QToolButton, text="Dupl.", tooltip="Duplicate current tab", change=self.duplicate_tab), Qt.TopRightCorner)
		self.tabs.tabCloseRequested.connect(self.close_tab)
		self.tabs.tabBarDoubleClicked.connect(self.renameoradd_tab)
		self.tabs.setCurrentIndex(main.config["series_currenttab"])
		self.tabs.currentChanged.connect(lambda x: main.config.__setitem__("series_currenttab", x))
		self.tabs.currentChanged.connect(lambda x: self.changed())

		self.set_state(main.config["series_references"])
		if not self.tabs.count():
			self.add_tab()
		
		main.signalclass.createdplots.connect(self.min_series)

	def close_tab(self, index):
		if self.tabs.count() == 1:
			return
		tab = self.tabs.widget(index)
		tab.deleteLater()
		self.tabs.removeTab(index)
	
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
			self.tabs.removeTab(i)
		for tabdata in values:
			self.add_tab(tabdata)
	
	def get_state(self):
		values = []
		for i in range(self.tabs.count()):
			tmp = self.tabs.widget(i).get_values()
			values.append(tmp)
		return(values)
	
	def changed(self):
		main.config["series_references"] = self.get_state()

class CatalogueWindow(EQDockWidget):
	def __init__(self, parent=None):
		super().__init__(parent)
		self.setWindowTitle("Catalogue of newly assigned Lines")
		
		mainwidget = QGroupBox()
		layout = QVBoxLayout()
		self.setWidget(mainwidget)
		mainwidget.setLayout(layout)
		
		button_savecatalogue = QQ(QToolButton, text="Save", change=lambda: main.save_lines_lin())
		button_delete = QQ(QToolButton, text="X", change=main.catalogue_table_delete)
		button_deleteall = QQ(QToolButton, text="Del All", change=lambda: main.catalogue_table_delete(all=True))
		button_addrow = QQ(QToolButton, text="+", change=lambda x: addemptyrow_inplace(main.new_df, self.catalogueModel))
		button_tight = QQ(QToolButton, text="Tight", change=lambda x: self.catalogueTable.resizeColumnsToContents())
		checkbox_appendonsave = QQ(QCheckBox, "flag_appendonsave", text="Append", tooltip="Append to file if checked or overwrite content if unchecked")
		errorlabel = QQ(QLabel, text="Default Uncertainty: ", tooltip="Positive Values are absolute Values, -1 for obs-calc, -2 for dialog, -3 for StdDev from Fit")
		errorfield = QQ(QDoubleSpinBox, "fit_error", range=(-3, PYQT_MAX), singlestep=main.config["fit_uncertaintystep"])

		headers = ["U1", "U2", "U3", "U4", "U5", "U6", "L1", "L2", "L3", "L4", "L5", "L6", "Freq", "Unc.", "Weight", "Comment"]
		str_columns = [headers.index("Comment")]

		self.catalogueTable = QTableView()
		self.catalogueModel = CustomTableModel(main.new_df, headers, ["filename"])
		self.catalogueTable.setModel(self.catalogueModel)
		self.catalogueTable.resizeColumnsToContents()
		
		main.signalclass.assignment.connect(self.scroll_bottom)
		
		main.config.register("series_qns", self.update_columns_visibility)
		self.update_columns_visibility()

		buttonsBox = QHBoxLayout()
		[buttonsBox.addWidget(button_delete), buttonsBox.addWidget(button_deleteall), buttonsBox.addWidget(button_addrow),
		 buttonsBox.addStretch(2), buttonsBox.addWidget(checkbox_appendonsave), buttonsBox.addWidget(button_savecatalogue),
		 buttonsBox.addWidget(button_tight)]
		layout.addLayout(buttonsBox)
		layout.addWidget(self.catalogueTable)
		buttonsBox = QHBoxLayout()
		[buttonsBox.addWidget(errorlabel), buttonsBox.addWidget(errorfield), buttonsBox.addStretch(2)]
		layout.addLayout(buttonsBox)
	
	def scroll_bottom(self):
		self.catalogueTable.selectRow(len(self.catalogueModel.data)-1)
		self.catalogueModel.update()
		self.catalogueTable.scrollToBottom()
	
	def update_columns_visibility(self):
		qns = main.config["series_qns"]
		for i in range(6):
			self.catalogueTable.setColumnHidden(i,   i>=qns)
			self.catalogueTable.setColumnHidden(i+6, i>=qns)

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
		tmp = self.log_area.toPlainText()
		tmp = tmp.split("\n")
		if len(tmp)-1 > main.config["flag_logmaxrows"]:
			self.log_area.setText("\n".join(tmp[-main.config["flag_logmaxrows"]:]))
		
		self.log_area.append(text)
		sb = self.log_area.verticalScrollBar()
		sb.setValue(sb.maximum())

class QuoteWindow(EQDockWidget):
	def __init__(self, parent=None):
		super().__init__(parent)
		self.setWindowTitle("Quote")
		
		mainwidget = QGroupBox()
		layout = QVBoxLayout()
		self.setWidget(mainwidget)
		mainwidget.setLayout(layout)
		
		self.quote = QQ(QLabel, wordwrap=True, align=Qt.AlignCenter)
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
		self.quote.setSizePolicy(QSizePolicy.Minimum, QSizePolicy.Minimum)




##
## Enhanced Window Class
##
class EQWidget(QWidget):
	def __init__(self, id, parent=None):
		self.id = id
		super().__init__(parent)

		geometry = main.config.get(f"windowgeometry_{self.id}")
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

class CreditsWindow(EQWidget):
	def __init__(self, id, parent=None):
		super().__init__(id, parent)
		self.setWindowTitle("Credits")

		global CREDITSSTRING
		layout = QVBoxLayout()
		layout.addWidget(QQ(QLabel, text=CREDITSSTRING, align=Qt.AlignCenter, wordwrap=True, minHeight=300, minWidth=500))
		self.setLayout(layout)

class FilesWindow(EQWidget):
	def __init__(self, id, actions, title, type, lock, parent=None):
		super().__init__(id, parent)
		self.setWindowTitle(title)
		self.actions = actions
		self.type = type
		self.lock = lock
		self.files_dict = {}
		main.signalclass.updatewindows.connect(self.update)
		self.layout = None
		self.scrollarea = None
		self.update()

	def gui_top(self):
		pass

	def gui_bottom(self):
		pass

	def gui_init(self):
		if self.layout == None:
			vbox = QVBoxLayout()

			self.scrollarea = QScrollArea()
			widget = QWidget()
			self.layout = QVBoxLayout()

			self.scrollarea.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
			self.scrollarea.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
			self.scrollarea.setWidgetResizable(True)

			widget.setLayout(self.layout)
			self.scrollarea.setWidget(widget)
			vbox.addWidget(self.scrollarea)

			self.setLayout(vbox)

			self.layout_files = QGridLayout()
			self.layout_actions = QHBoxLayout()
			self.layout.addLayout(self.layout_actions)
			self.layout_actions.addWidget(QQ(QPushButton, text="Load", change=self.add_file))
			self.layout_actions.addWidget(QQ(QPushButton, text="Add", change=lambda: self.add_file(True)))
			if "reread" in self.actions:
				self.layout_actions.addWidget(QQ(QPushButton, text="Reread All", change=lambda x: self.reread_file()))
			self.layout_actions.addWidget(QQ(QPushButton, text="Delete All", change=lambda x: self.delete_file()))
			self.layout_actions.addStretch(1)
			self.setMinimumHeight(200)
			self.setMinimumWidth(600)
			self.layout.addLayout(self.layout_files)
			self.layout.addStretch(1)

		self.row_id = 0
		self.gui_top()
		for file in main.config[f"files_{self.type}"]:
			self.add_row(file, self.actions)
		self.gui_bottom()

	def add_row(self, file, actions):
		layout = self.layout_files
		
		color = main.config[f"files_{self.type}"][file].get("color", None)
		hidden = main.config[f"files_{self.type}"][file].get("hidden")
		
		rowdict = {
			"label":			QQ(QLabel, text=file, enabled=not hidden),
			"colorinput":		QQ(QLineEdit, text=color, maxWidth=200, change=lambda x, file=file: self.change_color(file, inp=True)),
			"colorpicker":		QQ(QPushButton, text="CP", change=lambda x, file=file: self.change_color(file), stylesheet=f"background-color: {rgbt_to_trgb(color)}"),
			"invert":			QQ(QPushButton, text="Invert", change=lambda x, file=file: self.invert_file(file)),
			"hide":				QQ(QPushButton, text="Show" if hidden else "Hide", change=lambda x, file=file: self.hide_file(file)),
			"delete":			QQ(QPushButton, text="X", change=lambda x, file=file: self.delete_file(file)),
			"reread":			QQ(QPushButton, text="Reread", change=lambda x, file=file: self.reread_file(file)),
		}

		layout.addWidget(rowdict["label"], self.row_id, 0)

		for col_id, action in enumerate(self.actions):
			layout.addWidget(rowdict[action], self.row_id, col_id)
	
		self.files_dict[file] = rowdict
		self.row_id += 1

	def update(self):
		if self.scrollarea:
			tmp = (self.scrollarea.verticalScrollBar().value(), self.scrollarea.horizontalScrollBar().value())
		else:
			tmp = (0, 0)
		for key, value in self.files_dict.items():
			for widget in value.values():
				widget.deleteLater()
		self.files_dict = {}
		self.gui_init()
		self.scrollarea.verticalScrollBar().setValue(tmp[0])
		self.scrollarea.horizontalScrollBar().setValue(tmp[1])

	@synchronized_d(locks["axs"])
	def hide_file(self, file, inp=False):
		hidden = main.config[f"files_{self.type}"][file].get("hidden", False)
		hidden = not hidden
		main.config[f"files_{self.type}"][file]["hidden"] = hidden

		if hidden:
			self.files_dict[file]["label"].setEnabled(False)
			self.files_dict[file]["hide"].setText("Show")
		else:
			self.files_dict[file]["label"].setEnabled(True)
			self.files_dict[file]["hide"].setText("Hide")

		main.signalclass.updateplot.emit()

	@threading_d
	@working_d
	@synchronized_d(locks["exp_df"])
	@synchronized_d(locks["cat_df"])
	@synchronized_d(locks["lin_df"])
	def delete_file(self, file=None):
		df = main.return_df(self.type)
		
		with self.lock:
			if file == None:
				main.config[f"files_{self.type}"].clear()
				df.drop(df.index, inplace=True)
			else:
				if file in main.config[f"files_{self.type}"]:
					del main.config[f"files_{self.type}"][file]
				df.drop(df[df["filename"]==file].index, inplace=True)
		
		main.load_file(self.type, keep_old=True, do_QNs=False)

	@threading_d
	@working_d
	@synchronized_d(locks["exp_df"])
	@synchronized_d(locks["cat_df"])
	@synchronized_d(locks["lin_df"])
	def reread_file(self, file=None):
		if not file:
			for file in main.config[f"files_{self.type}"]:
				main.load_file(self.type, add_files=[file], keep_old=True)
		else:
			main.load_file(self.type, add_files=[file], keep_old=True)

class CatWindow(FilesWindow):
	def __init__(self, id, parent=None):
		type = "cat"
		title = "Cat Files"
		actions = ["label", "colorinput", "colorpicker", "hide", "delete", "reread"]
		lock = locks["cat_df"]
		super().__init__(id, actions, title, type, lock, parent)

	def add_file(self, keep_old=False):
		main.load_file("cat", add_files=True, keep_old=keep_old)

	def gui_top(self):
		for file, title, color in zip(("__starassigned__", "__current__"), ("Assigned Markers", "Current Transition"), (main.config["color_lin"], main.config["color_cur"])):
			rowdict = {	
				"label":			QQ(QLabel, text=title),
				"colorinput":		QQ(QLineEdit, text=color, maxWidth=200, change=lambda x, file=file: self.change_color(file, inp=True)),
				"colorpicker":		QQ(QPushButton, text="CP", change=lambda x, file=file: self.change_color(file), stylesheet=f"background-color: {rgbt_to_trgb(color)}"),
			}

			for col_id, action in enumerate(("label", "colorinput", "colorpicker")):
				self.layout_files.addWidget(rowdict[action], self.row_id, col_id)
	
			self.files_dict[file] = rowdict
			self.row_id += 1

	@synchronized_d(locks["axs"])
	def change_color(self, file, inp=False):
		if inp:
			color = self.files_dict[file]["colorinput"].text()
		else:
			color = QColorDialog.getColor(initial=QColor(rgbt_to_trgb(self.files_dict[file]["colorinput"].text())), options=QColorDialog.ShowAlphaChannel)
			if color.isValid():
				color = trgb_to_rgbt(color.name(QColor.HexArgb))
			else:
				return
		
		try:
			color = Color(color)
		except CustomError:
			return
			
		self.files_dict[file]["colorpicker"].setStyleSheet(f"background-color: {rgbt_to_trgb(color)}")
		if self.files_dict[file]["colorinput"].text() != color:
			self.files_dict[file]["colorinput"].setText(color)

		if file == "__starassigned__":
			main.config["color_lin"] = color
			for plot in main.plotwidget.axs["lin_plot"].flatten():
				plot.set_color(color)
		elif file == "__current__":
			main.config["color_cur"] = color
		else:
			main.config[f"files_{self.type}"][file]["color"] = color
		main.signalclass.updateplot.emit()

class ExpWindow(FilesWindow):
	def __init__(self, id, parent=None):
		type = "exp"
		title = "Exp Files"
		actions = ["label", "colorinput", "colorpicker", "reread", "hide", "invert", "delete"]
		lock = locks["exp_df"]
		super().__init__(id, actions, title, type, lock, parent)

	def add_file(self, keep_old=False):
		main.load_file("exp", add_files=True, keep_old=keep_old)

	def gui_top(self):
		color = main.config.get('color_exp')
		file = "__defaultexpcolor__"
		rowdict = {	
			"label":			QQ(QLabel, text="Initial Color"),
			"colorinput":		QQ(QLineEdit, text=color, maxWidth=200, change=lambda x, file=file: self.change_color(file, inp=True)),
			"colorpicker":		QQ(QPushButton, text="CP", change=lambda x, file=file: self.change_color(file), stylesheet=f"background-color: {rgbt_to_trgb(color)}"),
		}

		for col_id, action in enumerate(("label", "colorinput", "colorpicker")):
			self.layout_files.addWidget(rowdict[action], self.row_id, col_id)

		self.files_dict[file] = rowdict
		self.row_id += 1

	def change_color(self, file, inp=False):
		if inp:
			color = self.files_dict[file]["colorinput"].text()
		else:
			color = QColorDialog.getColor(initial=QColor(rgbt_to_trgb(self.files_dict[file]["colorinput"].text())), options=QColorDialog.ShowAlphaChannel)
			if color.isValid():
				color = trgb_to_rgbt(color.name(QColor.HexArgb))
			else:
				return
		
		try:
			color = Color(color)
		except CustomError:
			return
		
		self.files_dict[file]["colorpicker"].setStyleSheet(f"background-color: {rgbt_to_trgb(color)}")
		if self.files_dict[file]["colorinput"].text() != color:
			self.files_dict[file]["colorinput"].setText(color)

		if file == "__defaultexpcolor__":
			main.config["color_exp"] = color
		else:
			main.config[f"files_{self.type}"][file]["color"] = color
		main.signalclass.updateplot.emit()

	@synchronized_d(locks["exp_df"])
	def invert_file(self, file):
		main.config[f"files_{self.type}"][file]["invert"] = not main.config[f"files_{self.type}"][file].get("invert", False)
		tmp_df = main.exp_df
		tmp_df.loc[tmp_df.filename == file, "y"] = -tmp_df.loc[tmp_df.filename == file, "y"]
		main.signalclass.updateplot.emit()

class LinWindow(FilesWindow):
	def __init__(self, id, parent=None):
		type = "lin"
		title = "Lin Files"
		actions = ["label", "reread", "hide", "delete"]
		lock = locks["lin_df"]
		super().__init__(id, actions, title, type, lock, parent)

	def add_file(self, keep_old=False):
		main.load_file("lin", add_files=True, keep_old=keep_old)

	def gui_top(self):
		file = "__lin_own_df__"
		hidden = main.config["flag_hidecatalogue"]
		rowdict = {
			"label":			QQ(QLabel, text="Newly Assigned Lines"),
			"hide":				QQ(QPushButton, text="Show" if hidden else "Hide", change=lambda x, file=file: self.hide_new_df(file)),
		}

		for col_id, action in enumerate(("label", "hide")):
			self.layout_files.addWidget(rowdict[action], self.row_id, col_id)
	
		self.files_dict[file] = rowdict
		self.row_id += 1

	@synchronized_d(locks["axs"])
	def hide_new_df(self, file, inp = False):
		hidden = not main.config.get("flag_hidecatalogue", False)
		main.config["flag_hidecatalogue"] = hidden
		if hidden:
			self.files_dict[file]["button_hide"].setText("Show")
			self.files_dict[file]["label"].setDisabled(True)
		else:
			self.files_dict[file]["hide"].setText("Hide")
			self.files_dict[file]["label"].setDisabled(False)
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
		layout.addWidget(self.plotWidget)

		tmp_layout = QHBoxLayout()
		layout.addLayout(tmp_layout)

		tmp_layout.addWidget(QQ(QLabel, text="Function: "))
		tmp_layout.addWidget(QQ(QComboBox, "lineshapewindow_lineshape", items=["Gauss", "Lorentz", "Voigt"]))

		tmp_layout.addWidget(QQ(QLabel, text="Gauss: "))
		tmp_layout.addWidget(QQ(QDoubleSpinBox, "lineshapewindow_gauss"))
		
		tmp_layout.addWidget(QQ(QLabel, text="Lorentz: "))
		tmp_layout.addWidget(QQ(QDoubleSpinBox, "lineshapewindow_lorentz"))
		
		tmp_layout.addWidget(QQ(QLabel, text="Derivative: "))
		tmp_layout.addWidget(QQ(QSpinBox, "lineshapewindow_derivative", range=(0, 2)))
		
		tmp_layout.addWidget(QQ(QLabel, text="Scaling: "))
		tmp_layout.addWidget(QQ(QDoubleSpinBox, "lineshapewindow_scaling_factor", range=(-PYQT_MAX, PYQT_MAX)))
		tmp_layout.addWidget(QQ(QSpinBox, "lineshapewindow_scaling_exponent", prefix="*10^", range=(-PYQT_MAX, PYQT_MAX)))

		tmp_layout.addStretch()

		main.config.register(("lineshapewindow_lineshape", "lineshapewindow_gauss", "lineshapewindow_lorentz", "lineshapewindow_derivative", "lineshapewindow_scaling_factor", "lineshapewindow_scaling_exponent"), self.plotWidget.update_plot)

	def activateWindow(self):
		self.plotWidget.from_current_plot()
		super().activateWindow()

class BlendedLinesWindow(EQWidget):
	def __init__(self, id, parent=None):
		super().__init__(id, parent)
		self.setWindowTitle("Blended Lines Window")

		self.peaks = []
		self.cid = None

		class CustomPlotWidget(ProtPlot):
			def gui(self):
				super().gui()
				self.fit_line = self.ax.plot([], [], color = main.config["blendedlineswindow_color_total"])[0]
				self.cat_line = None
				self.plot_parts = []

			def fit_peaks(self):
				with locks["currThread"]:
					ownid = threading.current_thread().ident
				
				try:
					self.parent.label.setText("<span style='color:#eda711;'>Working ...</span>")
					peaks = self.parent.peaks.copy()
					profile = main.config["blendedlineswindow_lineshape"]
					derivative = main.config["blendedlineswindow_derivative"]
					polynomrank = main.config["blendedlineswindow_polynom"]+1
					fixedwidth = main.config["blendedlineswindow_fixedwidth"]
					now = 2 if profile == "Voigt" else 1
					noa = 2 + now * (not fixedwidth)

					fitfunction = lambda x, *args, fixedwidth=fixedwidth: self.parent.fitfunction(x, profile, derivative, polynomrank, *args, fixedwidth=fixedwidth)
					fitfunction_withoutbaseline = lambda x, *args, fixedwidth=fixedwidth: self.parent.fitfunction(x, profile, derivative, 0, *args, fixedwidth=fixedwidth)
					
					for part in self.plot_parts:
						part.remove()
					self.plot_parts = []

					xrange = [self.center-self.width/2, self.center+self.width/2]

					breakpoint(ownid, self.fit_thread_id)

					df_exp = main.get_visible_data("exp", xrange=xrange)
					exp_xs = df_exp["x"].to_numpy()
					exp_ys = df_exp["y"].to_numpy()

					breakpoint(ownid, self.fit_thread_id)

					if self.cat_line != None:
						self.cat_line.remove()
						self.cat_line = None
						
					df_cat = main.get_visible_data("cat", xrange=xrange)
					if len(df_cat) != 0:
						cat_xs = df_cat["x"].to_numpy()
						cat_ys = df_cat["y"].to_numpy()
						if len(exp_ys):
							cat_ys = cat_ys*np.max(exp_ys)/np.max(cat_ys)
						
						segs = (((cat_xs[i],0),(cat_xs[i],cat_ys[i])) for i in range(len(cat_xs)))
						colors = create_colors(df_cat, main.config["files_cat"])
						self.cat_line = self.ax.add_collection(matplotlib.collections.LineCollection(segs, colors=colors))

					breakpoint(ownid, self.fit_thread_id)

					xs = []
					ys = []
					ws = []
					if len(exp_ys) and (polynomrank + len(peaks)):
						p0 = []
						bounds = [[],[]]
						
						y0 = 4*(np.amax(exp_ys)-np.amin(exp_ys))
						w0 = xrange[1] - xrange[0]
						wmax = main.config["blendedlineswindow_maxfwhm"] or w0
						for peak in peaks:
							x, y, x_rel = peak

							if not xrange[0] < x < xrange[1]:
								x = self.center + x_rel
								if not xrange[0] < x < xrange[1]:
									x = sum(xrange)/2
								y = y0/2
							elif not 0 < y < y0:
								y = y0/2
							
							xs.append((x, *xrange))
							ys.append((y, 0, y0))
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
						
						try:
							popt, pcov = optimize.curve_fit(fitfunction, exp_xs, exp_ys, p0=p0, bounds=bounds)
						except Exception as E:
							popt, pcov = optimize.curve_fit(fitfunction, exp_xs, exp_ys, p0=p0, bounds=bounds)
						perr = np.sqrt(np.diag(pcov))
						res_xs = np.linspace(xrange[0], xrange[1], main.config["blendedlineswindow_xpoints"])
						res_ys = fitfunction(res_xs, *popt)
					else:
						popt = [0]*(noa*len(peaks)+polynomrank+2*fixedwidth)
						perr = [0]*(noa*len(peaks)+polynomrank+2*fixedwidth)
						res_xs = np.linspace(xrange[0], xrange[1], main.config["blendedlineswindow_xpoints"])
						res_ys = res_xs*0
					
					breakpoint(ownid, self.fit_thread_id)
					
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
						self.plot_parts.append(self.ax.plot(res_xs, tmp_ys, color=main.config["blendedlineswindow_color_total"], alpha=main.config["blendedlineswindow_transparency"])[0])
						
						opt_param.append( tmp_params+list(peaks[i]) )
						err_param.append( tmp_errors+list(peaks[i]) )
					
					self.parent.params = opt_param, err_param, profile, noa, now

					self.plot_parts.append(self.ax.scatter([x[0] for x in opt_param], [x[1] for x in opt_param], color=main.config["blendedlineswindow_color_points"]))

					if polynomrank > 0 and main.config["blendedlineswindow_showbaseline"]:
						self.plot_parts.append(self.ax.plot(res_xs, np.polyval(popt[-polynomrank:], res_xs-self.center), color=main.config["blendedlineswindow_color_baseline"])[0])

					breakpoint(ownid, self.fit_thread_id)

					main.signalclass.blwfit.emit()
					self.parent.label.setText("Ready")
					super().update_plot()
				except CustomError as E:
					pass
				except Exception as E:
					self.parent.label.setText("<span style='color:#ff0000;'>Fit failed</span>")
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

		row1 = (
		  QQ(QLabel, text="Function: "), 
		  QQ(QComboBox, "blendedlineswindow_lineshape", items=("Gauss", "Lorentz", "Voigt")),
		  QQ(QLabel, text="Derivative: "),
		  QQ(QSpinBox, "blendedlineswindow_derivative", range=(0, 2), minWidth=50),
		  QQ(QLabel, text="Max FWHM: "),
		  QQ(QDoubleSpinBox, "blendedlineswindow_maxfwhm", range=(0, PYQT_MAX)),
		  QQ(QLabel, text="Transparency: "),
		  QQ(QDoubleSpinBox, "blendedlineswindow_transparency", range=(0, 1), minWidth=50, singlestep=0.1),
		  None,
		)
	
		row2 = (
		  QQ(QLabel, text="Baseline Rank: "),
		  QQ(QSpinBox, "blendedlineswindow_polynom", range=(-1, PYQT_MAX)),
		  QQ(QCheckBox, "blendedlineswindow_showbaseline", text="Show Baseline"),
		  QQ(QCheckBox, "blendedlineswindow_fixedwidth", text="All Same Width"),
		  None,
		)
		
		self.label = QQ(QLabel, text="Ready", textFormat=Qt.RichText)
		row3 = (
		  QQ(QPushButton, text="Del All", change=lambda x: self.del_peak(-1)),
		  QQ(QPushButton, text="Update", change=lambda x: self.plotWidget.update_plot()),
		  self.label,
		  None,
		)

		for row in (row1, row2, row3):
			tmp_layout = QHBoxLayout()
			layout.addLayout(tmp_layout)
			for widget in row:
				if widget is None:
					tmp_layout.addStretch()
				else:
					tmp_layout.addWidget(widget)

		self.table = QTableWidget()
		self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
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
		opt_param, err_param, function, noa, now = self.params
		opt_param.sort(key=lambda x: x[0])
		table = self.table
		table.setRowCount(0)
		table.setColumnCount(7)
		table.setHorizontalHeaderLabels(["Action", "Frequency", "Amplitude", "FWHM Gauss", "FWHM Lorentz", "Other QNs", "Delete"])
		for params, errparam in zip(opt_param, err_param):
			currRowCount = table.rowCount()
			table.insertRow(currRowCount)
			table.setCellWidget(currRowCount, 0, QQ(QPushButton, text="Assign", change=lambda x, xpos=params[0], error=errparam[1]: self.pre_assign(xpos, error)))
			table.setItem(currRowCount, 1, QTableWidgetItem(f'{params[0]:.4f}'))
			table.setItem(currRowCount, 2, QTableWidgetItem(f'{params[1]:.4f}'))
			table.setItem(currRowCount, 3, QTableWidgetItem(f'{params[2] if function != "Lorentz" else  0:.4f}'))
			table.setItem(currRowCount, 4, QTableWidgetItem(f'{params[1+now] if function != "Gauss" else  0:.4f}'))
			table.setCellWidget(currRowCount, 5, QQ(QPushButton, text="Assign other QNs", change=lambda x, xpos=params[0], error=errparam[1]: self.pre_assign(xpos, error, oqns=True)))
			table.setCellWidget(currRowCount, 6, QQ(QPushButton, text="Delete", change=lambda x, ind=currRowCount: self.del_peak(i=ind)))
		table.resizeColumnsToContents()

	def pre_assign(self, x, error, oqns=False):
		index = self.plotWidget.i
		
		dict_ = {
			"x": x,
			"error": error,
			"xpre": 0,
		}
		
		for i in range(6):
			dict_[f"qnu{i+1}"] = np.iinfo(np.int64).min
			dict_[f"qnl{i+1}"] = np.iinfo(np.int64).min
		
		if oqns:
			visible_cat_data = main.get_visible_data("cat")
			dialog = QNsDialog(x, visible_cat_data)
			dialog.exec_()

			if dialog.result() == 1:
				dict_.update(dialog.save())
			else:
				return
		else:
			reference = main.config["series_references"][index[1]]
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
		self.outputTable.setEditTriggers(QAbstractItemView.NoEditTriggers)

		vertLayout = QHBoxLayout()
		leftLayout = QGridLayout()
		rightLayout = QVBoxLayout()
		layout.addWidget(QQ(QLabel, wordwrap=True, text="Series Finder allows to find the strongest (unassigned) predicted transitions."))

		rightLayout.addWidget(QQ(QLabel, text="Allowed Transitions"))
		rightLayout.addWidget(QQ(QCheckBox, "seriesfinderwindow_atype", text="a-type"))
		rightLayout.addWidget(QQ(QCheckBox, "seriesfinderwindow_btype", text="b-type"))
		rightLayout.addWidget(QQ(QCheckBox, "seriesfinderwindow_ctype", text="c-type"))
		rightLayout.addStretch(1)

		leftLayout.addWidget(QQ(QLabel, text="Frequency Range:"), 0, 0)
		leftLayout.addWidget(QQ(QLabel, text="Start Frequency: "), 1, 0)
		leftLayout.addWidget(QQ(QLineEdit, "seriesfinderwindow_start"), 1, 1)
		leftLayout.addWidget(QQ(QLabel, text="Stop Frequency: "), 2, 0, 1, 1)
		leftLayout.addWidget(QQ(QLineEdit, "seriesfinderwindow_stop"), 2, 1, 1, 1)
		leftLayout.addWidget(QQ(QLabel, text="Number of Results: "), 3, 0, 1, 1)
		leftLayout.addWidget(QQ(QSpinBox, "seriesfinderwindow_results", range=(0, PYQT_MAX)), 3, 1, 1, 1)
		leftLayout.addWidget(QQ(QLabel, text="Additional Condition: "), 4, 0, 1, 1)
		leftLayout.addWidget(QQ(QLineEdit, "seriesfinderwindow_condition"), 4, 1, 1, 1)
		leftLayout.addWidget(QQ(QCheckBox, "seriesfinderwindow_onlyunassigned", text = "Only unassigned Lines"), 5, 1, 1, 1)
		leftLayout.addWidget(QQ(QPushButton, text="Run", change=lambda x: self.run()), 6, 1, 1, 1)

		vertLayout.addLayout(leftLayout)
		vertLayout.addLayout(rightLayout)
		vertLayout.addStretch(1)
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
			tmp_condition.append(f"(abs(qnu2-qnl2) == 0 and abs(qnu3-qnl3) == 1)")
		if main.config["seriesfinderwindow_btype"]:
			tmp_condition.append(f"(abs(qnu2-qnl2) == 1 and abs(qnu3-qnl3) == 1)")
		if main.config["seriesfinderwindow_ctype"]:
			tmp_condition.append(f"(abs(qnu2-qnl2) == 1 and abs(qnu3-qnl3) == 0)")
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

		tmp_cat_df["y"] = np.log(tmp_cat_df["y"])/np.log(10)
		tmp_cat_df = tmp_cat_df.nlargest(nor, "y")

		if tmp_min and tmp_max:
			xrange = f"in the range from {tmp_min} to {tmp_max}"
		else:
			xrange = "in the total range"
		message = f"The {nor} most intense predicted transitions {unassigned} {xrange} are shown below."
		self.messageLabel.setText(message)
		self.messageLabel.setHidden(False)

		table = self.outputTable
		headers = ["Start", "Intensity", "Frequency"] + qns_visible

		table.setRowCount(0)
		table.setColumnCount(len(headers))
		table.setHorizontalHeaderLabels(headers)

		for index, row in tmp_cat_df.iterrows():
			currRowCount = table.rowCount()
			table.insertRow(currRowCount)
			table.setCellWidget(currRowCount,0, QQ(QPushButton, text="Start", change=lambda x, crow=row: self.startHere(crow)))
			table.setItem(currRowCount, 1, QTableWidgetItem(f'{row["y"]:.4f}'))
			table.setItem(currRowCount, 2, QTableWidgetItem(f'{row["x"]:.4f}'))

			for i, column in enumerate(qns_visible):
				tmp = row[column]
				if tmp == np.iinfo(np.int64).min:
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
					text = f"{entry[val]:.4f}".rstrip("0").rstrip(".")
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

		self.noq = main.config["series_qns"]
		self.table = QTableWidget()
		self.cols = ["x", "y", "dist"] + [f"qn{ul}{i+1}" for ul in ("u", "l") for i in range(6)] + ["filename"]
		self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
		self.table.setRowCount(0)
		self.table.setColumnCount(len(self.cols)+1)
		self.table.setHorizontalHeaderLabels(["Y/N", "x", "y", "Dist"] +  [f"{ul}{i+1}" for ul in ("U", "L") for i in range(6)] + ["Filename"])
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
				"comment":	"",
				"filename":	"__lin_own_df__",
			}

			for i in range(6):
				tmp_dict[f"qnu{i+1}"] = np.iinfo(np.int64).min
				tmp_dict[f"qnl{i+1}"] = np.iinfo(np.int64).min
			
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
		
		layout.addWidget(QQ(QPlainTextEdit, value=command, change=lambda: self.update_pipe_command(), placeholder="Write your command line commando here"))
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
		self.peaks_df = pd.DataFrame(columns=["x", "y", "go to"])
		main.signalclass.updatetable.connect(self.update_table)

		class CustomPlotWidget(ProtPlot):
			def gui(self):
				super().gui()
				self.peaks_line = self.ax.scatter([], [], color = main.config["peakfinderwindow_peakcolor"], marker="*")

			def update_plot(self):
				xmin, xmax = self.center-self.width/2, self.center+self.width/2
				tuples = list(filter(lambda x: xmin < x[0] < xmax, self.parent.peaks))
				tuples = tuples if len(tuples)!=0 else [[None,None]]
				self.peaks_line.set_offsets(tuples)
				super().update_plot()

		layout = QVBoxLayout()
		self.setLayout(layout)

		self.plotWidget = CustomPlotWidget(parent=self)
		layout.addWidget(self.plotWidget)

		tmp_layout = QHBoxLayout()
		layout.addLayout(tmp_layout)

		self.run_button = QQ(QPushButton, text="Run Peakfinder", change=lambda x: self.find_peaks())
		tmp_layout.addWidget(self.run_button)
		tmp_layout.addWidget(QQ(QPushButton, text="Export Peaks", change=lambda x: self.export_peaks()))

		tmp_layout.addWidget(QQ(QLabel, text="Uncertainty: "))
		uncert_input = QQ(QDoubleSpinBox, "peakfinderwindow_width", range=(0, PYQT_MAX))
		tmp_layout.addWidget(QQ(QCheckBox, "peakfinderwindow_onlyunassigned", text="Only unassigned Lines", change=uncert_input.setEnabled))
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
			input_min  = QQ(QDoubleSpinBox, range=(0, PYQT_MAX))
			input_max  = QQ(QDoubleSpinBox, range=(0, PYQT_MAX))

			tmp_layout.addWidget(label, i+1, 0)
			tmp_layout.addWidget(input_min, i+1, 1)
			if key == "distance":
				items = (("Min", "Off"))
				input_min.setRange(1, PYQT_MAX)
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
						input_min.setValue(value[0])
						input_type.setCurrentIndex(0)
				except:
					pass
			else:
				if key == "distance":
					input_type.setCurrentIndex(1)
				else:
					input_type.setCurrentIndex(2)
			self.change_type(key)

		self.infolabel = QQ(QLabel, wordwrap=True, text="Press 'Run Peakfinder' to find peaks.")
		layout.addWidget(self.infolabel)

		self.table = QTableView()
		self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
		self.table.setHidden(True)
		self.table.clicked.connect(self.go_to)
		self.tableModel = CustomTableModel(self.peaks_df, ["x", "y", "Go to"], editable=False)
		self.table.setModel(self.tableModel)
		layout.addWidget(self.table)
		
		main.signalclass.peakfinderend.connect(self.update_table)

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
		self.run_button.setEnabled(False)
		self.infolabel.setText("Finding peaks ...")

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
			xs = exp_df["x"].to_numpy()
			ys = exp_df["y"].to_numpy()

		peaks = signal.find_peaks(ys, **kwargs)[0]

		if len(peaks) == 0:
			tuples = []
		else:
			px, py = xs[peaks], ys[peaks]
			tuples = list(zip(px, py))

		self.peaks = tuples
		self.plotWidget.update_plot()

		peak_df = self.get_peaks()
		if len(peak_df) > main.config["peakfinderwindow_maxentries"]:
			self.infolabel.setText(f"Found {len(peak_df)} peaks. Only the highest {main.config['peakfinderwindow_maxentries']} entries are shown in the table for a better performance. You can see all peaks by exporting the peaks or increasing the maximum number of displayed peaks in the configuration.")
			peak_df.sort_values("y", ascending=False, inplace=True)
			peak_df = peak_df.head(main.config["peakfinderwindow_maxentries"])
		else:
			self.infolabel.setText(f"Found {len(peak_df)} peaks.")

		self.peaks_df.drop(self.peaks_df.index, inplace=True)
		self.peaks_df.reset_index(drop=True, inplace=True)
		peak_df.reset_index(drop=True, inplace=True)
		for i in range(len(peak_df)):
			self.peaks_df.loc[i] = peak_df.loc[i]

		main.signalclass.peakfinderend.emit()

	@synchronized_d(locks["peaks"])
	def update_table(self):
		self.table.setHidden(False)
		self.tableModel.update()
		self.run_button.setEnabled(True)

	@synchronized_d(locks["peaks"])
	def export_peaks(self):
		peak_df = self.get_peaks()[["x", "y"]]
		fname = QFileDialog.getSaveFileName(None, 'Choose file to save peaks to',"","CSV Files (*.csv);;All Files (*)")[0]
		if fname:
			peak_df.to_csv(fname, header=None, index=None, sep='\t')

	@synchronized_d(locks["peaks"])
	def get_peaks(self):
		peaks = self.peaks
		peak_df = pd.DataFrame(peaks, columns=["x", "y"])

		if main.config["peakfinderwindow_onlyunassigned"]:
			u = main.config["peakfinderwindow_width"]
			peak_df["assigned"] = False

			assigned_xs = main.get_visible_data("lin")["x"]
			for x in assigned_xs:
				peak_df.loc[((x-u) < peak_df["x"]) & (peak_df["x"] < (x+u)), "assigned"] = True

			peak_df.drop(peak_df[peak_df["assigned"]==True].index, inplace=True)
			peak_df.sort_values("y", ascending=False, inplace=True)
			del peak_df["assigned"]

		peak_df.sort_values("y", ascending=False, inplace=True)
		peak_df.reset_index(drop=True, inplace=True)
		peak_df["go to"] = "Go here"

		return(peak_df)

	@synchronized_d(locks["peaks"])
	def go_to(self, event):
		for idx in self.table.selectionModel().selectedIndexes():
			row_number = idx.row()
			column_number = idx.column()

		if column_number == 2 and type(row_number) == int and len(self.peaks_df) > row_number:
			x = self.peaks_df.loc[row_number, "x"]
			self.plotWidget.center = x
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

		self.catalogueTable = QTableView()
		self.catalogueModel = CustomTableModel(main.ser_df, list(main.ser_df.columns))
		self.catalogueTable.setModel(self.catalogueModel)
		self.catalogueTable.resizeColumnsToContents()

		class CustomSeriesSelector(SeriesSelector):
			def changed(self):
				super().changed()
				main.config["seriesfitwindow_series"] = self.state
		
		self.series_selector = CustomSeriesSelector(self, main.config["seriesfitwindow_series"])
		self.output_area = QQ(QPlainTextEdit, readonly=True, minHeight=50, placeholder="This is the output field and is read-only.")

		layout = QVBoxLayout()
		buttonsBox = QHBoxLayout()
		buttonsBox.addWidget(QQ(QToolButton, text="X", change=lambda x: self.catalogueTableDelete()))
		buttonsBox.addWidget(QQ(QToolButton, text="+", change=lambda x: addemptyrow_inplace(main.ser_df, self.catalogueModel)))
		buttonsBox.addWidget(QQ(QToolButton, text="Copy Lines From Main Window", change=lambda x: self.copyMWLines()))
		buttonsBox.addWidget(QQ(QCheckBox, "seriesfitwindow_greedy", text="Add Newly Assigned Lines here"))
		buttonsBox.addStretch(1)
		buttonsBox.addWidget(QQ(QToolButton, text="Tight", change=lambda x: self.catalogueTable.resizeColumnsToContents()))
		layout.addLayout(buttonsBox)
		layout.addWidget(self.catalogueTable)
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
		self.catalogueTable.selectRow(len(self.catalogueModel.data)-1)
		self.catalogueModel.update()
		self.catalogueTable.scrollToBottom()

	@synchronized_d(locks["ser_df"])
	def catalogueTableDelete(self):
		selected = [index.row() for index in self.catalogueTable.selectionModel().selectedRows()]
		for index in sorted(selected, reverse = True):
			main.ser_df.drop(index, inplace =True)
			main.ser_df.reset_index(inplace = True, drop = True)
		self.catalogueModel.update()

	def copyMWLines(self):
		main.ser_df.reset_index(drop=True, inplace=True)
		main.new_df.reset_index(drop=True, inplace=True)
		length = len(main.ser_df)
		tmp_df = main.new_df[main.ser_df.columns]
		
		for i in range(len(main.new_df)):
			tmp_dict = tmp_df.loc[i]
			tmp_dict = {key: main.ser_dtypes[key](value) for key, value in tmp_dict.items()}
			main.ser_df.loc[length+i] = tmp_dict
		self.catalogueModel.update()

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
		df.sort_values("x", inplace = True)
		
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
		
		self.pred_qns = [np.concatenate((qnu+incrs*i, qnl+incrs*i)) for i in range(main.config["seriesfitwindow_maxprediction"])]
		self.pred_xs  = [self.function(qns, *popt) for qns in self.pred_qns]

		tmp = "\n".join([f"{name} : {value:.4f}" for name, value in zip(self.fitparams, popt)])
		self.writelog(f"Succeeded, parameters were determined as \n{tmp}")

	def show_pred_freqs(self):
		reference = main.config["series_references"][main.config["series_currenttab"]]
		reference["method"] = "List"
		reference["list"] = {"qns": self.pred_qns, "xs": self.pred_xs, "i0": 0} # @Luis: Check if this format works for the qns
		main.config["series_qns"] = self.noq
		
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
		layout.addWidget(QQ(QPlainTextEdit, "residualswindow_query", maxHeight=40, placeholder="Query text to filter shown lines. Use qnu1, ..., qnu6 and qnl1, ..., qnl6 for the quantum numbers."))
		layout.addWidget(QQ(QPlainTextEdit, "residualswindow_colorinput", maxHeight=40, placeholder="Enter custom color and query to color specific lines differently. E.g. enter '#ff0000; qnu1 < 20' to color all transitions with the first upper quantum number below 20 red."))
		
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

		return(df)
	
	def plot_residuals(self):
		self.update_button.setDisabled(True)
		main.app.processEvents()
		try:
			message = []
			
			df = self.get_residuals()
			df["obs_calc"] = df["x_lin"] - df["x_cat"]
			df["color"] = main.config["residualswindow_defaultcolor"]
			message.append(f"Found {len(df)} entries matching your query.")

			colorquerytext = main.config["residualswindow_colorinput"].split("\n")
			for row in colorquerytext:
				if row.strip():
					color, query = row.split(";")
					indices = df.query(query).index
					df.iloc[indices, df.columns.get_loc("color")] = color
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
				self.ax.set_xlim([np.min(xs), np.max(xs)])
				y_range = [np.min(ys), np.max(ys)]
				self.ax.set_ylim(y_range[0]-main.config["plot_ymargin"]*(y_range[1]-y_range[0]), y_range[1]+main.config["plot_ymargin"]*(y_range[1]-y_range[0]))
			self.fig.canvas.draw_idle()
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
			self.fig.canvas.draw_idle()
			
			if click and noq:
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
		
		self.fig = figure.Figure(dpi=main.config["plot_dpi"])
		self.plot_canvas = FigureCanvas(self.fig)
		
		self.ax = self.fig.subplots()
		self.ax.ticklabel_format(useOffset=False)
		
		self.points = self.ax.scatter([], [], color=[], marker=".")
		self.annot = self.ax.annotate("", xy=(0,0), xytext=(5,5), textcoords="offset points", color="black", ha="center", va="bottom", bbox=dict(boxstyle="round", fc="w"))
		self.annot.set_visible(False)
		self.fig.canvas.mpl_connect("motion_notify_event", self.on_hover)
		
		self.mpl_toolbar = NavigationToolbar2QT(self.plot_canvas, self)
		
		self.fname = None
		self.df = None
		
		layout = QVBoxLayout()
		self.setLayout(layout)
		layout.addWidget(self.plot_canvas, 6)
		layout.addWidget(self.mpl_toolbar)
		hlayout = QHBoxLayout()
		hlayout.addWidget(QQ(QPushButton, text="Open", change=self.load_file))
		self.file_label = QQ(QLabel, text="No File loaded")
		hlayout.addWidget(self.file_label)
		hlayout.addWidget(QLabel("x-axis: "))
		hlayout.addWidget(QQ(QLineEdit, "residualswindow_xvariable", placeholder="Choose the x-variable, e.g. qn1"))
		hlayout.addWidget(QLabel("y-axis: "))
		hlayout.addWidget(QQ(QLineEdit, "residualswindow_yvariable", placeholder="Choose the x-variable, e.g. egy"))
		hlayout.addWidget(QQ(QCheckBox, "residualswindow_autoscale", text="Autoscale on Update"))
		layout.addLayout(hlayout)
		layout.addWidget(QQ(QPlainTextEdit, "energylevelswindow_query", maxHeight=40, placeholder="Query text to filter shown levels. Use qn1, ..., qn6 for the quantum numbers."))
		layout.addWidget(QQ(QPlainTextEdit, "energylevelswindow_colorinput", maxHeight=40, placeholder="Enter custom color and query to color specific lines differently. E.g. enter '#ff0000; qn1 < 20' to color all levels with the first quantum number below 20 red."))

		buttonslayout = QHBoxLayout()
		layout.addLayout(buttonslayout)
		buttonslayout.addStretch(1)
		self.update_button = QQ(QPushButton, text="Update", change=self.plot_energylevels)
		buttonslayout.addWidget(self.update_button)
		buttonslayout.addStretch(1)
	
	def get_energylevels(self, fname):
		## Egy Dataframe Format
		_egy_df_columns = ['iblk', 'indx', 'egy', 'err', 'pmix', 'we', ':', 'qn1', 'qn2', 'qn3', 'qn4', 'qn5', 'qn6']
		_egy_df_dtypes = [np.int64, np.int64, np.float64, np.float64, np.float64, np.int64, str]+[np.int64]*6
		_egy_df_widths = [0,6,11,29,47,58,63,64,67,70,73,76,79,82]
		
		dtypes_dict = {_egy_df_columns[i]: _egy_df_dtypes[i] for i in range(len(_egy_df_columns))}
		
		data = []
		with open(fname, "r") as file:
			for line in file:
				if line.strip() == "" or line.startswith("#"):
					continue
				data.append([column_to_numeric(line[i:j]) for i,j in zip(_egy_df_widths[:-1], _egy_df_widths[1:])])
		data = pd.DataFrame(data)
		data.columns = _egy_df_columns
		data = data.astype(dtypes_dict)
		data.reset_index(drop=True, inplace=True)
		return(data)
	
	def load_file(self):
		fname = QFileDialog.getOpenFileName(None, 'Choose Egy File to load',"")[0]
		if fname:
			self.fname = fname
			self.plot_energylevels()
			self.file_label.setText(os.path.split(fname)[1])
	
	def plot_energylevels(self):
		if self.fname == None:
			return
		self.update_button.setDisabled(True)
		main.app.processEvents()
		try:
			df = self.get_energylevels(self.fname)
			query = main.config["energylevelswindow_query"]
			if query:
				df = df.query(query, inplace=True)
			df["color"] = main.config["energylevelswindow_defaultcolor"]

			colorquerytext = main.config["energylevelswindow_colorinput"].split("\n")
			for row in colorquerytext:
				if row.strip():
					color, query = row.split(";")
					df.iloc[df.query(query).index, df.columns.get_loc("color")] = color
			
			self.df = df
			
			yvariable = main.config["energylevelswindow_xvariable"].strip() or "egy"
			xvariable = main.config["energylevelswindow_xvariable"].strip() or "qn1"
			ys = df.eval(yvariable).to_numpy()
			xs = df.eval(xvariable).to_numpy()
			colors = df["color"].to_numpy()
			tuples = list(zip(xs,ys))
			tuples = tuples if len(tuples)!=0 else [[None,None]]
			self.points.set_offsets(tuples)
			self.points.set_color(colors)
			if len(xs) and self.autoscale.isChecked():
				self.ax.set_xlim([np.min(xs), np.max(xs)])
				y_range = [np.min(ys), np.max(ys)]
				self.ax.set_ylim(y_range[0]-main.config["plot_ymargin"]*(y_range[1]-y_range[0]), y_range[1]+main.config["plot_ymargin"]*(y_range[1]-y_range[0]))
			self.fig.canvas.draw_idle()
		except:
			main.notification("There was an error in your Energy Levels window inputs")
			raise
		finally:
			self.update_button.setDisabled(False)
	
	def on_hover(self, event):
		if event.inaxes == self.ax and isinstance(self.df, pd.DataFrame):
			cont, ind = self.points.contains(event)
			if cont:
				self.annot.xy = self.points.get_offsets()[ind["ind"][0]]
				tmp_levels = self.df.iloc[ind["ind"]]
				text = []
				for i, row in tmp_levels.iterrows():
					qns = row[self.qns_visible].to_numpy()
					text.append(",".join(str(int(row[f"qn{i+1}"])) for i in range(main.config["series_qns"])))
				text = "\n".join(text)
				self.annot.set_text(text)
				self.annot.set_visible(True)
			else:
				self.annot.set_visible(False)
			self.fig.canvas.draw_idle()

class SpectraResolverWindow(EQWidget):
	def __init__(self, id, parent=None):
		super().__init__(id, parent)
		self.setWindowTitle("Spectra Resolver")
		self.list_ = QListWidget()
		self.list_.setDragDropMode(QAbstractItemView.InternalMove)
		
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
		main.signalclass.updatewindows.connect(self.fill_list)
	
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
		self.solveoverlap_core(fname_savefile)
	
	@threading_d
	def solveoverlap_core(self, fname_savefile):
		self.apply_order_button.setEnabled(False)
		main.app.processEvents()
		self.label.setText("Working")
		try:
			with locks["exp_df"]:
				df = main.exp_df.copy()
			df["keep"] = 1
			i_keep = df.columns.get_loc("keep")
			
			ranked_files = [self.list_.item(i).text() for i in range(self.list_.count())]
			results = []

			for i, fname in enumerate(ranked_files):
				self.label.setText(f"Working on file {i+1} of {len(ranked_files)}")
				min, max = main.config["files_exp"][fname]["xrange"]
				i_start, i_stop = df["x"].searchsorted(min, side="left"), df["x"].searchsorted(max, side="right")
				
				tmp_df = df.iloc[i_start:i_stop]
				df.iloc[i_start:i_stop, i_keep] = 0
				tmp_df = tmp_df[tmp_df["filename"] == fname]
				results.append(tmp_df.copy())
				df = df[df["keep"] == 1]
				
			self.label.setText(f"Working on saving the results")
			df = pd.concat(results, ignore_index=True)
			df.sort_values("x", inplace=True, kind="merge")			
			df[["x", "y"]].to_csv(fname_savefile, header=None, index=None, sep=chr(main.config["flag_separator"]))
		finally:
			self.label.setText(f"Ready")
			self.apply_order_button.setEnabled(True)

class ReportWindow(EQWidget):
	def __init__(self, id, parent=None):
		super().__init__(id, parent)
		self.setWindowTitle("Report Window")
		
		layout = QVBoxLayout()
		self.setLayout(layout)
		
		layout.addWidget(QQ(QPlainTextEdit, "reportwindow_query", maxHeight=40, placeholder="Query text to filter assignments. Use qnu1, ..., qnu6 and qnl1, ..., qnl6 for the quantum numbers."))
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
		tmplayout2.addWidget(QQ(QDoubleSpinBox, "figurewindow_width", range=(0, PYQT_MAX)))
		tmplayout2.addWidget(QQ(QLabel, text=" x "))
		tmplayout2.addWidget(QQ(QDoubleSpinBox, "figurewindow_height", range=(0, PYQT_MAX)))
		tmplayout2.addWidget(QQ(QComboBox, "figurewindow_unit", options=("cm", "inch")))
		tmplayout.addLayout(tmplayout2, 0, 1)
		
		row_index += 1
		tmplayout.addWidget(QQ(QLabel, text="DPI: "), row_index, 0)
		tmplayout.addWidget(QQ(QDoubleSpinBox, "figurewindow_dpi", range=(0, PYQT_MAX)), row_index, 1)
		
		row_index += 1
		tmplayout.addWidget(QQ(QLabel, text="Font-Size: "), row_index, 0)
		tmplayout.addWidget(QQ(QDoubleSpinBox, "figurewindow_fontsize", range=(0, PYQT_MAX)), row_index, 1)
		
		row_index += 1
		tmplayout.addWidget(QQ(QLabel, text="x-Label: "), row_index, 0)
		tmplayout.addWidget(QQ(QLineEdit, "figurewindow_xlabel"), row_index, 1)
		
		row_index += 1
		tmplayout.addWidget(QQ(QLabel, text="x-Label Padding: "), row_index, 0)
		tmplayout.addWidget(QQ(QDoubleSpinBox, "figurewindow_xlabelpadding", range=(-PYQT_MAX, PYQT_MAX)), row_index, 1)
		
		row_index += 1
		tmplayout.addWidget(QQ(QLabel, text="y-Label: "), row_index, 0)
		tmplayout.addWidget(QQ(QLineEdit, "figurewindow_ylabel"), row_index, 1)
		
		row_index += 1
		tmplayout.addWidget(QQ(QLabel, text="y-Label Padding: "), row_index, 0)
		tmplayout.addWidget(QQ(QDoubleSpinBox, "figurewindow_ylabelpadding", range=(-PYQT_MAX, PYQT_MAX)), row_index, 1)
		
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
			QKeySequence(QKeySequence.ZoomIn),
			self.view,
			activated=self.zoom_in,
		)
		
		QShortcut(
			QKeySequence(QKeySequence.ZoomIn),
			self.view,
			activated=self.zoom_in,
		)

		QShortcut(
			QKeySequence(QKeySequence.ZoomOut),
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
			
			xpos, allqns = main.plotwidget.get_positions(return_qns=True)
			widths = main.plotwidget.get_widths()
			
			xmin = np.zeros(xpos.shape)
			xmax = np.zeros(xpos.shape)
			
			for i in range(rows):
				for j in range(cols):
					ax = self.axs[i, j]
					ax.yaxis.set_visible(False)
					ax.xaxis.set_visible(False)
					offset = main.plotwidget.get_offset((i, j))
					x, width = xpos[i, j]+offset, widths[i, j]/2
					xrange = (x-width, x+width)
					ax.set_xlim(*xrange)
					
					xmin[i, j], xmax[i, j] = xrange
			
			# set ticks for bottom row
			i = -1
			for j in range(cols):
				ax = self.axs[i, j]
				ax.xaxis.set_visible(True)
				ticks = np.linspace(xmin[i, j], xmax[i, j], main.config["plot_ticks"])
				if main.config["plot_offsetticks"] == 1:
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
			
			# set data and set y-ranges of data
			scaling = main.config["plot_yscale"]
			
			dataframes = [main.get_visible_data(x) for x in ("exp", "cat", "lin")]
			files_dicts = [main.config[x] for x in ("files_exp", "files_cat", "files_lin")]
			ann_dict = main.config["plot_annotations_dict"]

			for i in range(rows):
				for j in range(cols):
					ax = self.axs[i, j]

					qns = allqns[i, j]
					color = matplotlib.rcParams['text.color']
					
					if main.config["figurewindow_annotation"] == "Custom":
						reference = main.config["series_references"][j]
						if reference["method"] == "Transition":
							qn_dict = {f"qn{ul}{i+1}": qns[0 + (ul=="l")][i] for ul in "ul" for i in range(len(qns[0]))}
						else:
							qn_dict = {}
						text = main.config["figurewindow_customannotation"].format(**{"x": xpos[i, j], "i": i, "j": j, "rows": rows, "cols": cols, **qn_dict})
					elif main.config["series_annotate_xs"] or qns is None:
						text = f"{{:{main.config['series_annotate_fmt']}}}".format(xpos[i, j])
					else:
						text = f"{', '.join([str(qn) if qn != np.iinfo(np.int64).min else '-' for qn in qns[0]])} ← {', '.join([str(qn) if qn != np.iinfo(np.int64).min else '-' for qn in qns[1]])}"
						
						lin_df = main.get_visible_data("lin")
						if len(lin_df.query(" and ".join([f"qnu{i+1} == {qn}" for i, qn in enumerate(qns[0])] + [f"qnl{i+1} == {qn}" for i, qn in enumerate(qns[1])]))):
							color = main.config["color_lin"]
				
					ax.text(**ann_dict, s=text, transform=ax.transAxes, color=color)
					
					
					for datatype, dataframe, files in zip(("exp", "cat", "lin"), dataframes, files_dicts):
						minindex = dataframe["x"].searchsorted(xmin[i, j], side="left")
						maxindex = dataframe["x"].searchsorted(xmax[i, j], side="right")
						
						dataframe = dataframe.iloc[minindex:maxindex].copy()

						xs = dataframe["x"].to_numpy()
						ys = dataframe["y"].to_numpy() if datatype != "lin" else 0*xs
						
						if datatype == "exp":
							if scaling == "Per Plot":
								if len(dataframe):
									yrange_exp = [dataframe["y"].min(), dataframe["y"].max()]
								else:
									yrange_exp = [-1, 1]
							ax.plot(xs, ys, **main.config["figurewindow_expkwargs"])
						elif datatype == "cat":
							if scaling == "Per Plot":
								if len(dataframe):
									yrange_cat = [dataframe["y"].min(), dataframe["y"].max()]
								else:
									yrange_cat = [-1, 1]
								ys = ys*yrange_exp[1]/yrange_cat[1]
							elif scaling in ["Global", "Custom"]:
								ys = ys*main.config["plot_expcat_factor"]*10**main.config["plot_expcat_exponent"]
							segs = (((xs[i], 0),(xs[i], ys[i])) for i in range(len(xs)))
							colors = create_colors(dataframe, files, xpos[i, j])
							coll = matplotlib.collections.LineCollection(segs, **{"colors": colors, **main.config["figurewindow_catkwargs"]})
							ax.add_collection(coll)
						elif datatype == "lin":
							ax.scatter(xs, ys, **{"color": main.config["color_lin"], **main.config["figurewindow_linkwargs"]})

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

		scrollarea.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
		scrollarea.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
		scrollarea.setWidgetResizable(True)

		tmp_layout = QHBoxLayout()
		tmp_layout.addWidget(QQ(QPushButton, text="Save as default", change=lambda: main.saveoptions()))
		completer = QCompleter(main.config.keys())
		completer.setCaseSensitivity(Qt.CaseInsensitive)
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
		if not converter:
			return
		else:
			converter = converter[1]
		input, oklab, label = self.widgets[key]

		try:
			if converter in (dict, list, tuple):
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




##
## Dialogs
##
class QNsDialog(QDialog):
	def __init__(self, frequency, df):
		super().__init__()
		QShortcut("Esc", self).activated.connect(lambda: self.predone(0))

		qns = [f"qnu{i+1}" for i in range(6)]+[f"qnl{i+1}" for i in range(6)]
		self.res = {key: np.iinfo(np.int64).min for key in qns}

		self.setWindowTitle(f"Choose QNs for transition at {frequency}")
		self.resize(main.config["qnsdialog_width"], main.config["qnsdialog_height"])
		self.frequency = frequency

		layout = QVBoxLayout()
		self.setLayout(layout)
		layout.addWidget(QQ(QLabel, wordwrap=True, text="Enter the quantum numbers in the input fields or choose one of the transitions from the table."))
		layout.addSpacing(10)

		noq = main.config["series_qns"]
		self.sbs = {}
		qnslayout = QGridLayout()
		for i in range(noq):
			widget = QQ(QSpinBox, value=0, range=(0, PYQT_MAX))
			qnslayout.addWidget(widget, 1, i)
			self.sbs[f"qnu{i+1}"] = widget
		
		for i in range(noq):
			widget = QQ(QSpinBox, value=0, range=(0, PYQT_MAX))
			qnslayout.addWidget(widget, 2, i)
			self.sbs[f"qnl{i+1}"] = widget
		
		for i in range(noq):
			lab = QLabel(f"QN {i+1}")
			qnslayout.addWidget(QQ(QLabel, text=f"QN{i+1}"), 0, i)

		qnslayout.setColumnStretch(7, 1)
		layout.addLayout(qnslayout)

		cols = ["dist", "x"] + qns
		table = QTableWidget()
		table.setColumnCount(len(cols)+1)
		table.setHorizontalHeaderLabels(["Assign"] + cols)
		table.setRowCount(0)
		table.setEditTriggers(QAbstractItemView.NoEditTriggers)
		
		tmp_df = df.copy()
		tmp_df["dist"] = frequency-tmp_df["x"]
		tmp_df["absdist"] = abs(tmp_df["dist"])
		tmp_df.sort_values(by=["absdist"], inplace=True)
		tmp_df.reset_index(drop=True, inplace=True)

		for i, row in tmp_df.head(100).iterrows():
			currRowCount = table.rowCount()
			table.insertRow(currRowCount)
			for j, col in enumerate(cols):
				val = f'{row[col]:.4f}'.rstrip("0").rstrip(".")
				table.setItem(currRowCount, j+1, QTableWidgetItem(val))
			tmpd = {key: row[key] for key in qns}
			tmpd["xpre"] = row["x"]
			table.setCellWidget(currRowCount, 0, QQ(QPushButton, text="Assign", change=lambda x, tmpd=tmpd: self.table_save(tmpd)))
		
		for i in range(6):
			table.setColumnHidden(i+3, i>=noq)
			table.setColumnHidden(i+9, i>=noq)
		
		table.resizeColumnsToContents()
		layout.addWidget(table)

		buttons = QDialogButtonBox.Ok | QDialogButtonBox.Cancel
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
		cursor.movePosition(QTextCursor.End)
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
			main.config["pipe_commands"].append(["New Tab", "", True, True, True])
		elif self.tabs.widget(index) != 0:
			text, ok = QInputDialog().getText(self, "Tab Name","Enter the Tabs Name:")
			if ok and text:
				self.tabs.setTabText(index, text)
				main.config["pipe_commands"][index][0] = text

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
		})
		
		if initial_values:
			self.state.update(initial_values)
		
		plotwidgets = {}
		self.methods = ("Transition", "List", "Expression")
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
		width_blend = QQ(QDoubleSpinBox, "series_blendwidth", range=(0, PYQT_MAX))
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
		self.xsTable.setEditTriggers(QAbstractItemView.NoEditTriggers)
		self.xsTable.setHorizontalHeaderLabels(["#", "Frequency"])
		layout.addWidget(self.xsTable)
		if self.state["list"]["xs"]:
			self.load_xs_list(values = self.state["list"]["xs"])
		
		button_open = QQ(QToolButton, text="Open List", change=lambda x: self.load_xs_list(temp=False))
		button_write = QQ(QToolButton, text="Write List", change=lambda x: self.load_xs_list(temp=True))

		label_startat = QLabel("Start at Index: ")
		spinbox_startat = QQ(QSpinBox, value=self.state["list"]["i0"], range=(0, PYQT_MAX), singlestep=1, change=lambda: self.state["list"].__setitem__("i0", spinbox_startat.value()))
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
		self.input_N0 = QQ(QSpinBox, value=self.state["expression"]["N0"], range=(0, PYQT_MAX), change=lambda: self.state["expression"].__setitem__("N0", self.input_N0.value()))

		layout.addWidget(self.input_expr)
		hbox = QHBoxLayout()
		[hbox.addWidget(label_N0), hbox.addWidget(self.input_N0), hbox.addStretch(1), hbox.addWidget(button_apply)]
		layout.addLayout(hbox)

		plotwidgets["Expression"].setLayout(layout)

	def get_values(self):
		return(self.state)
	
	def changed(self):
		pass
		
	def update(self, newstate):
		self.state.update(newstate)
		self.series_selector.state = self.state["transition"]
		self.series_selector.set_state()

	def load_xs_list(self, values=None, temp=False):
		if values:
			self.state["list"]["xs"] = values
			xs = values
		
		elif temp:
			line, ok = QInputDialog().getMultiLineText(self, "Specify Custom List",
			"Write list here (delimiters are all whitespace characters, comma and semicolon):")
			if not ok or not line:
				return
			
			xs = []
			tmp_xs = re.split('; |, |\s', line)
			for x in tmp_xs:
				try:
					xs.append(float(x))
				except ValueError:
					main.notification(f"<span style='color:#eda711;'>WARNING</span>: Could not convert the string '{x}' to a numerical value.")
			
			self.state["list"] = {
				"qns":		None,
				"xs":		xs,
				"i0":	0,
			}
		
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
						
						tmp = re.split('; |, |\s', line)
						for x in tmp:
							try:
								xs.append(float(x))
							except ValueError:
								main.notification(f"<span style='color:#eda711;'>WARNING</span>: Could not convert the string '{x}' to a numerical value.")
			
			self.state["list"] = {
				"qns":		None,
				"xs":		xs,
				"i0":	0,
			}
		
		

		table = self.xsTable
		table.setRowCount(0)
		i=0
		for i, x in enumerate(xs):
			currRowCount = table.rowCount()
			table.insertRow(currRowCount)
			table.setItem(currRowCount, 0, QTableWidgetItem(f"{i}"))
			table.setItem(currRowCount, 1, QTableWidgetItem(f"{x:.4f}"))
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
		self.qnus    = [QQ(QSpinBox, maxWidth=40, range=(0, PYQT_MAX), singlestep=1, change=lambda x: self.changed()) for x in range(6)]
		self.qnls    = [QQ(QSpinBox, maxWidth=40, range=(0, PYQT_MAX), singlestep=1, change=lambda x: self.changed()) for x in range(6)]
		self.qnincrs = [QQ(QCheckBox, text="Incr", maxWidth=40, change=lambda x: self.changed()) for x in range(6)]
		self.qndiffs = [QQ(QSpinBox, maxWidth=40, range=(-PYQT_MAX, PYQT_MAX), singlestep=1, change=lambda x: self.changed()) for x in range(6)]

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

		layout.addWidget(self.togglediff, 4, 6, 1, 1)
		
		layout.addWidget(self.incqns, 1, 6, 1, 2)
		layout.addWidget(self.decqns, 2, 6, 1, 2)

		layout.addWidget(self.incnoq, 0, 6)
		layout.addWidget(self.decnoq, 0, 7)

		layout.setRowStretch(6, 10)
		layout.setColumnStretch(8, 10)
		
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
class ProtPlot(QWidget):
	def __init__(self, parent=None, i=None, onrange=False):
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
			self.span_selector = matplotlib.widgets.SpanSelector(self.ax, lambda vmax, vmin: self.on_range(vmax, vmin, self.i), 'horizontal')

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

		self.exp_line = self.ax.plot([], [], color=main.config["color_exp"])[0]

	def from_current_plot(self, update=True):
		self.i = main.plotwidget.get_current_plot()
		xrange = main.plotwidget.axs["ax"][self.i].get_xlim()
		self.center = sum(xrange)/2
		self.width = xrange[1] - xrange[0]

		if update:
			self.update_plot()

	@synchronized_d(locks["exp_df"])
	def update_plot(self):
		exp_df = main.get_visible_data("exp", xrange=(self.center-self.width/2, self.center+self.width/2), binning=True)

		self.exp_xs = exp_df["x"].to_numpy()
		self.exp_ys = exp_df["y"].to_numpy()

		self.exp_line.set_data(self.exp_xs, self.exp_ys)
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
		}
		
		
		for keys, function in shortcuts_dict.items():
			QShortcut(keys, self.parent).activated.connect(function)
	
	@synchronized_d(locks["axs"])
	def on_range(self, xmin, xmax, index):
		axrange = self.axs["ax"][index].get_xlim()
		if xmax == xmin or xmax > axrange[1] or xmin < axrange[0]:
			return
		shift = (QApplication.keyboardModifiers() == Qt.ShiftModifier)
		if self.onrange == "Assign" and (shift or main.config["fit_alwaysfit"]):
			xmiddle, xuncert = main.plotwidget.fit_peak(xmin, xmax, index)
			xpre = main.plotwidget.self.cache_positions[index]
			
			dict_ = {
				"x":		xmiddle,
				"error":	xuncert,
				"xpre":		xpre
			}
			
			reference = main.config["series_references"][index[1]]
			if reference["method"] == "Transition":
				qns = main.plotwidget.get_qns(reference["transition"])[index[0]]
				for i, qnu, qnl in zip(range(6), qns[0], qns[1]):
					dict_[f"qnu{i+1}"] = qnu
					dict_[f"qnl{i+1}"] = qnl
	
			if not main.config["fit_alwaysassign"]:
				main.plotwidget.lastassignment = dict_
			elif main.plotwidget.check_blends(index, dict_):
				return
			else:
				main.plotwidget.assign(index, dict_)
		else:
			self.center = (xmin+xmax)/2
			self.width = xmax-xmin
			self.update_plot()

class NotificationsBox(QWidget):
	def __init__(self):
		super().__init__()
		self.bg_color = QColor("#a5aab3")
		self.messages = []
		self.setWindowFlags(
			Qt.Window | Qt.Tool | Qt.FramelessWindowHint | 
			Qt.WindowStaysOnTopHint | Qt.X11BypassWindowManagerHint)

		self.setAttribute(Qt.WA_NoSystemBackground, True)
		self.setAttribute(Qt.WA_TranslucentBackground, True)

		self.setMinimumHeight(80)
		self.setMinimumWidth(300)
		self.setMaximumWidth(300)
		
		self.layout = QVBoxLayout()
		self.setLayout(self.layout)

		self.setStyleSheet("""
			color: white;
			background-color: #bf29292a;
		""")

		self._desktop = QApplication.instance().desktop()
		startPos = QPoint(self._desktop.screenGeometry().width() - self.width() - 10, 10)
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
		self.callbacks = pd.DataFrame(columns=["id", "key", "widget", "function"], dtype="object").astype({"id": np.uint})
		
	def __setitem__(self, key, value, widget=None):
		super().__setitem__(key, value)
		self.signal()
		
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
			self.callbacks = self.callbacks.append({
				"id":		0,
				"key":		key,
				"function":	function
			}, ignore_index=True)
	
	def register_widget(self, key, widget, function):
		ids = set(self.callbacks["id"])
		id = 1
		while id in ids:
			id += 1
		self.callbacks = self.callbacks.append({
			"id":		id,
			"key":		key,
			"widget":	widget,
			"function":	function
		}, ignore_index=True)
		widget.destroyed.connect(lambda x, id=id: self.unregister_widget(id))

	def unregister_widget(self, id):
		self.callbacks.drop(self.callbacks[self.callbacks["id"] == id].index, inplace=True)

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
	updateconfig  = pyqtSignal()
	updatewindows = pyqtSignal()
	updatetable   = pyqtSignal()
	assignment    = pyqtSignal()
	updateplot    = pyqtSignal()
	createdplots  = pyqtSignal()
	blwfit        = pyqtSignal()
	peakfinderend = pyqtSignal()
	writelog      = pyqtSignal(str)
	notification  = pyqtSignal(str)
	def __init__(self):
		super().__init__()

class QSpinBox(QSpinBox):
	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		self.setStepType(QAbstractSpinBox.AdaptiveDecimalStepType)
	
	def setSingleStep(self, value):
		self.setStepType(QAbstractSpinBox.DefaultStepType)
		super().setSingleStep(value)
	
	def setValue(self, value):
		if value < -PYQT_MAX or value > PYQT_MAX:
			value = 0
		return super().setValue(value)
		

class QDoubleSpinBox(QDoubleSpinBox):
	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		self.setDecimals(20)
		self.setStepType(QAbstractSpinBox.AdaptiveDecimalStepType)
	
	def setSingleStep(self, value):
		self.setStepType(QAbstractSpinBox.DefaultStepType)
		super().setSingleStep(value)
	
	def textFromValue(self, value):
		return(f"{value:.10f}".rstrip("0").rstrip("."))

	def valueFromText(self, text):
		return(np.float64(text))

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
		if role == Qt.DisplayRole or role == Qt.EditRole:
			value = self.data.iloc[index.row(), index.column()]
			dtype = self.data[self.data.columns[index.column()]].dtypes
			
			if isinstance(value, str):
				return(value)
			elif isinstance(value, (np.integer, int)):
				if value == np.iinfo(np.int64).min:
					return("")
				else:
					return(f"{{:{main.config['flag_tableformatint']}}}".format(value))
			elif np.isnan(value):
				return("")
			else:
				return(f"{{:{main.config['flag_tableformatfloat']}}}".format(value))

	def rowCount(self, index):
		return(self.data.shape[0])

	def columnCount(self, index):
		return(self.data.shape[1]-len(self.hiddencolumns))

	def headerData(self, section, orientation, role):
		if role == Qt.DisplayRole:
			if orientation == Qt.Horizontal:
				return str(self.headers[section])

			if orientation == Qt.Vertical:
				return str(self.data.index[section])

	def flags(self, index):
		if self.editable:
			return(Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable)
		else:
			return(Qt.ItemIsEnabled | Qt.ItemIsSelectable)

	def update(self):
		self.layoutChanged.emit()

	def setData(self, index, value, role):
		if not index.isValid():
			return False
		if role != Qt.EditRole:
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
					value = np.iinfo(np.int64).min
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

def csv_copypaste(self, event):
	if event.key() == Qt.Key_C and (event.modifiers() & Qt.ControlModifier):
		cells = sorted(self.selectedIndexes())
		output = []
		i = 0

		while i < len(cells):
			tmp = []
			row = cells[i].row()
			while i < len(cells) and cells[i].row() == row:
				tmp.append(cells[i].data())
				i += 1
			output.append("\t".join(tmp))
		output = "\n".join(output)
		QApplication.clipboard().setText(output)
	
	elif event.key() == Qt.Key_V and (event.modifiers() & Qt.ControlModifier):
		if QAbstractItemView.NoEditTriggers == self.editTriggers():
			return
		cells = sorted(self.selectedIndexes())
		if not cells:
			return
		text = QApplication.clipboard().text()
		data = [row.split("\t") for row in text.split("\n")]
		i_0, j_0 = cells[0].row(), cells[0].column()
		
		j_hidden = 0
		for i, row in enumerate(data):
			for j, value in enumerate(row):
				while self.isColumnHidden(j_0+j+j_hidden):
					j_hidden += 1
				self.model().setData(self.model().index(i_0+i, j_0+j+j_hidden), value, Qt.EditRole)
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

def create_colors(dataframe, files={}, xpos=None):
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
			return(np.iinfo(np.int64).min)
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
		dec_a = len(f"{a:.4f}".rstrip("0").split(".")[1])
		dec_o = len(f"{o:.4f}".rstrip("0").split(".")[1])
		if dec_a == dec_o:
			tick_labels.append(f"{a:.4f}".rstrip("0").rstrip("."))
		else:
			min_dec = max(dec_a, dec_o)
			tick_labels.append(f"{a:.4f}"[:-(4-min_dec)])
	return(tick_labels)

def except_hook(cls, exception, traceback):
	sys.__excepthook__(cls, exception, traceback)
	with open(f"{APP_TAG}.err", "a+", encoding="utf-8") as file:
		time_str = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
		file.write(f"{time_str}: \n{exception}\n{''.join(tb.format_tb(traceback))}\n\n")
	try:
		if not main.new_df.empty:
			main.save_lines_lin(f"{APP_TAG}.lin", force_noappend=True, force_lin=True)
		if main.config["flag_debug"]:
			main.self.notification(f"{exception}\n{''.join(tb.format_tb(traceback))}")
	except Exception as E:
		pass

def breakpoint(ownid, lastid):
	if ownid != lastid:
		raise CustomError()

def commandline(showdialog=True):
	if showdialog:
		dialog = ConsoleDialog()
		dialog.exec_()
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
	quotes = json.loads(quotes_str)
	quote = quotes[random.randint(0,len(quotes)-1)]
	quote = quote.replace("\n", "%0d%0a")
	webbrowser.open(f"mailto:bonah@ph1.uni-koeln.de?subject={APP_TAG} &body=%0d%0a%0d%0a%0d%0a%0d%0a{quote}")

def restart():
	dir = os.path.dirname(os.path.realpath(__file__))
	fname = os.path.join(dir, f"{APP_TAG}.llwp")
	main.saveproject(fname)
	main.saveoptions()
	if not main.new_df.empty:
		main.save_lines_lin(f"{APP_TAG}.lin", force_noappend=True, force_lin=True)
	os.execv(sys.executable, [sys.executable, sys.argv[0], fname])

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
		ys = np.gradient(ys)
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
			newvalues.append(np.iinfo(np.int64).min)
		else:
			newvalues.append("")
	df.loc[len(df.index)] = newvalues
	if model:
		model.update()

def exp_to_df(fname, sep="\t", xcolumn=0, ycolumn=1, invert=False, sort=True):
	data = pd.read_csv(fname, sep=sep, dtype=np.float64, header=None, engine="c", comment="#")
	column_names = [i for i in range(len(data.columns))]
	
	column_names[xcolumn if xcolumn in column_names else 0] = "x"
	column_names[ycolumn if ycolumn in column_names else 1] = "y"
	
	data.columns = column_names
	data = data[["x", "y",]]
	data["filename"] = fname
	
	if invert:
		data["y"] = -data["y"]
	return(data)

def shared_labels(fig, xlabel, ylabel, xlabelpad=15, ylabelpad=0, **kwargs):
	ax = fig.add_subplot(111, frameon=False)
	for side in ("top", "bottom", "left", "right"):
		ax.spines[side].set_color('none')
	ax.tick_params(axis="both", labelcolor='#fff0', top=False, bottom=False, left=False, right=False, zorder=-100)
	ax.set_xlabel(xlabel, labelpad=xlabelpad, **kwargs)
	ax.set_ylabel(ylabel, labelpad=ylabelpad, **kwargs)
	return(ax)

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
	"layout_theme":							["light", str],
	"layout_owntheme":						[{}, dict],
	"layout_limitannotationtosingleline":	[True, bool],
	"layout_mpltoolbar":					[False, bool],
	
	"color_exp":							["#000000", Color],
	"color_cat":							["#ff1a29", Color],
	"color_lin":							["#ff38fc", Color],
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
		
	"plot_dpi":								[100, float],
	"plot_annotations_dict":				[{"x": 1, "y": 1, "horizontalalignment": "right", "verticalalignment": "top"}, dict],
	"plot_font_dict":						[{"size":10}, dict],
	"plot_matplotlibkwargs":				[{"hspace": 0, "wspace": 0}, dict],
	"plot_coupled":							[True, bool],
	"plot_ymargin":							[0.1, float],
	"plot_annotate":						[True, bool],
	"plot_hover":							[True, bool],
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
	
	"series_qns":							[3, int],
	"series_annotate_xs":					[False, bool],
	"series_annotate_fmt":					[".4f", str],
	"series_blendwidth":					[1, float],
	"series_blenddialog":					[True, bool],
	"series_currenttab":					[1, int],
	"series_references":					[[], list],
		
	"flag_automatic_draw":					[True, bool],
	"flag_appendonsave":					[True, bool],
	"flag_hidecatalogue":					[False, bool],
	"flag_xcolumn":							[0, int],
	"flag_ycolumn":							[1, int],
	"flag_separator":						[9, int],
	"flag_debug":							[False, bool],
	"flag_alwaysshowlog":					[True, bool],
	"flag_extensions":						[{"exp": [".csv"], "cat": [".cat"], "lin": [".lin"], "project": [".llwp"]}, dict],
	"flag_autosetqns":						[True, bool],
	"flag_predictionformats":				[{}, dict, True],
	"flag_assignmentformats":				[{}, dict, True],
	"flag_assignmentsavefmt":				[{}, dict, True],
	"flag_loadfilesthreaded":				[True, bool],
	"flag_shownotification":				[True, bool],
	"flag_notificationtime":				[2000, int],
	"flag_showmainplotcontrols":			[True, bool],
	"flag_showmainplotwidth":				[True, bool],
	"flag_showmainplotrowscols":			[True, bool],
	"flag_referencenumberlocked":			[True, bool],
	"flag_logmaxrows":						[10000, int],
	"flag_tableformatint":					[".0f", str],
	"flag_tableformatfloat":				[".2f", str],
		
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
	
	"files_exp":							[{}, dict],
	"files_cat":							[{}, dict],
	"files_lin":							[{}, dict],
}

# Generate with following command, quotes is str of quotes.txt
# json.dumps(quotes.split("\n\n"))
quotes_str = '["Erstens kommt es anders und zweitens als man denkt\\n-Ekrem Bora", "In der Theorie willst du was ver\\u00e4ndern, aber in der Praxis geht alles schief.\\nIn der Theorie wollen wir zu viel, geben aber viel zu wenig, das ist kein fairer Deal.\\n-Anis M. Y. Ferchichi & Jochen Burchard", "Und bist du unten, dr\\u00fccken sie dich noch ein St\\u00fcck tiefer\\nNoch ein St\\u00fcck tiefer\\nNoch ein St\\u00fcck tiefer\\nNoch ein St\\u00fcck tiefer\\n-Anis M. Y. Ferchichi", "Da ting goes skrrrrahh, pap, pap, ka-ka-ka\\nSkidiki-pap-pap, and a pu-pu-pudrrrr-boom\\nSkya, du-du-ku-ku-dun-dun\\nPoom, poom, you dun know\\n-Michael Dapaah", "It wasn\'t me\\n-Orville Richard Burrell", "We are going so fast\\nBut time so slow\\n-Theory of Relativity", "D\\u00f6p, d\\u00f6p, d\\u00f6p, d\\u00f6p, d\\u00f6d\\u00f6d\\u00f6d\\u00f6d\\u00f6p\\nD\\u00f6, d\\u00f6d\\u00f6d\\u00f6p, d\\u00f6p, d\\u00f6d\\u00f6d\\u00f6d\\u00f6d\\u00f6p, d\\u00f6p\\n-Hans Peter Geerdes", "Skibadee, skibadanger\\nI am the rearranger\\n-Hans Peter Geerdes", "Respect to the Man in the Ice Cream Van\\n-Hans Peter Geerdes", "Hyper Hyper\\n-Hans Peter Geerdes", "If we could only slow the time\\nWe would have forever every night\\n-Don Pepijn Schipper", "Meine Stra\\u00dfenpoesie l\\u00f6st die Chaostheorie\\n-Mousa Amouei", "Die Parabel sie steigt, und zwar exponentiell\\n-Mohamed El Moussaoui", "Chuba chuba chuba chuba chuba chuba chubby.\\nI don\'t have any lines to go right here, so chuby Teletubby\\n-Marshall Bruce Mathers III", "Two things are infinite: The universe and human stupidity;\\nand I\\u2018m not sure about the universe\\n-Albert E", "Physics is like sex: sure, it may give some practical results, but that\'s not why we do it.\\n-Richard P. Feynman", "I do not think you can name many great inventions that have been made by married men.\\n-Nikola Tesla", "Those who are not shocked when they first come across quantum theory cannot possibly have understood it.\\n-Niels Bohr", "We\\u2019re not free in what we do, because we\\u2019re not free in what we want.\\n-Jonas Kahnwald", "What we know is a drop. What we don\\u2019t know is an ocean.\\n-Isaac Newton", "Das ist der Sound f\\u00fcr die echten M\\u00e4nner, die das hier h\\u00f6ren, wenn sie Pressluft h\\u00e4mmern\\n-Tarek Ebene, Nico Seyfrid & Maxim Dr\\u00fcner", "I accept that\\n-Chuck Marstein", "Many of life\'s failures are people who did not realize how close they were to success when they gave up\\n-Thomas A. Edison", "I find that the harder I work, the more luck I seem to have\\n-Thomas Jefferson", "Whether you think you can or you think you can\'t, you\'re right\\n-Henry Ford", "Life is never fair, and perhaps it is a good thing for most of us that it is not\\n-Oscar Wilde", "Only a life lived for others is a life worthwhile\\n-Albert Einstein", "Imagination is more important than knowledge\\n-Albert Einstein", "My mama always said, \\u2018Life was like a box of chocolates. You never know what you\\u2019re gonna get\\n-Forrest", "Before you marry a person, you should first make them use a computer with slow Internet to see who they really are\\n-Will Ferrell", "I don\\u2019t believe in astrology; I\\u2019m a Sagittarius and we\\u2019re skeptical\\n-Arthur C. Clarke", "If you think you are too small to make a difference, try sleeping with a mosquito\\n-Dalai Lama", "People who think they know everything are a great annoyance to those of us who do\\n-Isaac Asimov"]'

if __name__ == '__main__':
	main = Main()
	main.gui()