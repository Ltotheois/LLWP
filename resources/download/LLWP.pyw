#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Luis Bonah
# Description : Loomis-Wood Plot software

credits_string = """Author: Luis Bonah

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

app_tag = "LLWP"


import sys, os, re
import time, wrapt
import random, json
import threading, queue
import configparser
import traceback as tb
from io import StringIO
import numpy as np
import pandas as pd
from scipy import optimize
from scipy import special
import subprocess
import webbrowser

from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *

import matplotlib
from matplotlib import style, figure
from matplotlib.backends.backend_qt5agg import FigureCanvas, NavigationToolbar2QT

import warnings
warnings.simplefilter('ignore', np.RankWarning)

##
## Global Decorators
##
def stopwatch(func):
	def timed(*args, **kwargs):
		start_time = time.time()
		result = func(*args, **kwargs)
		stop_time = time.time()
		print(f"Executing {func.__name__} took {stop_time-start_time}.")
		return result
	return timed

def threading_d(func):
	def run(*args, **kwargs):
		if kwargs.get("add_files") != False and kwargs.get("reread") != True and type(kwargs.get("add_files")) != list:
			kwargs["add_files"] = QFileDialog.getOpenFileNames(None, 'Open File',)[0]
			if len(kwargs["add_files"])==0:
				return
		t = threading.Thread(target=func, args=args, kwargs=kwargs)
		t.start()
		return t
	return run

def synchronized(lock):
	@wrapt.decorator
	def _wrapper(wrapped, instance, args, kwargs):
		with lock:
			return wrapped(*args, **kwargs)
	return _wrapper

locks = {
	"exp_df" :					threading.RLock(),
	"lin_df" :					threading.RLock(),
	"ass_df" :					threading.RLock(),
	"cat_df" :					threading.RLock(),
	"ser_df" :					threading.RLock(),
	"windows" :					threading.RLock(),
	"axs" :						threading.RLock(),
	"pipe" :					threading.RLock(),
	"currThread" :				threading.RLock(),
}


##
## Main Window
##
class MainWindow(QMainWindow):
	def __init__(self, parent=None):
		super(MainWindow, self).__init__(parent)
		self.setFocusPolicy(Qt.StrongFocus)
		self.setWindowTitle(app_tag)
		self.setAcceptDrops(True)
		
		## Set Icon
		try:
			app.setWindowIcon(QIcon(f"{app_tag}.ico"))
			import ctypes
			ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(app_tag)
		except Exception as E:
			pass

		self.window_class = {
			"linesDialog_window" :		LinesWindow,
			"spectraDialog_window" :	SpectraWindow,
			"assignedDialog_window" :	AssignedWindow,
			"lineshape_window" :		LineshapeWindow,
			"blendedlines_window" :		BlendedlinesWindow,
			"pathfinder_window" :		PathfinderWindow,
			"pipe_window" :				PipeWindow,
			"series_fit_window" :		SeriesFitWindow,
			"credits_window" :			CreditsWindow,
		}
		self.open_windows = {key: None for key in self.window_class}
		self.load_options()
		self.main_widget.update_plots()
	
	
	##
	## Working Decorator
	##
	def working_d(func):
		def wrapper(self, *args, **kwargs):
			q = self.main_widget.working
			ind = self.main_widget.indicator
			if q.empty():
				ind.setText("Working...")
			q.put(1)
			try:
				res = func(self, *args, **kwargs)
			except Exception as E:
				exc_type, exc_value, exc_traceback = sys.exc_info()
				self.main_widget.notification(f"An error occurred: {str(E)}\n {''.join(tb.format_tb(exc_traceback))}")
				res = str(E)
			finally:
				q.get()
				q.task_done()
			if q.empty():
				ind.setText("Ready")
			return(res)
		return(wrapper)


	##
	## Top Level Functions
	##
	def save_options(self):
		self.config["plot_yrange_exp"]			= [self.main_widget.radio_yexp.buttons()[x].isChecked() for x in range(len(self.main_widget.radio_yexp.buttons()))].index(True)
		self.config["plot_yrange_lin"]			= [self.main_widget.radio_ylin.buttons()[x].isChecked() for x in range(len(self.main_widget.radio_ylin.buttons()))].index(True)
		self.config["plot_show_xaxis"]			= self.main_widget.checkbox_xaxis.isChecked()
		self.config["plot_coupled"]				= self.main_widget.checkbox_coupled.isChecked()
		self.config["plot_show_yaxis_exp"]		= self.main_widget.checkbox_show_yaxis_exp.isChecked()
		self.config["plot_show_yaxis_lin"]		= self.main_widget.checkbox_show_yaxis_lin.isChecked()
		self.config["plot_annotate"]			= self.main_widget.checkbox_annotate.isChecked()
		self.config["plot_hover_cutoff"]		= self.main_widget.checkbox_hover.isChecked()
		self.config["plot_number_of_plots"]		= len(self.main_widget._axs)
		self.config["plot_expyzoom"]			= self.main_widget.expScaleZoomInput.value()
		self.config["plot_linyzoom"]			= self.main_widget.linScaleZoomInput.value()

		self.config["series_number_of_plots"]	= self.main_widget.series_nop_input.text()
		self.config["series_width"]				= self.main_widget.series_width_input.text()
		self.config["series_ignore"]			= self.main_widget.series_ignore_checkbox.isChecked()
		self.config["series_J0"]				= self.main_widget.input_J0.text()
		self.config["series_expression"]		= self.main_widget.inputExpr.text()
		self.config["series_method"]			= self.main_widget._series_method
		self.config["series_transition"]		= json.dumps([	[x.value() for x in self.main_widget.QNUs],
																[x.value() for x in self.main_widget.QNLs],
																[x.isChecked() for x in self.main_widget.QNcbs],
															])

		self.config["fit_error"]				= self.main_widget.error_field.text()
		self.config["fit_alwaysfit"]			= self.always_fit_action.isChecked()
		self.config["fit_alwaysassign"]			= self.always_assign_action.isChecked()
		self.config["flag_append"]				= self.save_lines_append_action.isChecked()
		self.config["flag_automatic_draw"]		= self.automatic_draw_action.isChecked()
		self.config["flag_binning"]				= self.main_widget.binning_field.value()
		self.config["flag_alwaysshowlog"]		= self.always_show_log_action.isChecked()

		self.config["layout_options_boxes"]		= (not self.main_widget.topLeftGroupBox.isHidden()) and (not self.main_widget.topRightGroupBox.isHidden())
		self.config["layout_lines_box"]			= not self.main_widget.catalogueBox.isHidden()
		self.config["layout_log_box"]			= not self.main_widget.logBox.isHidden()
		self.config["layout_maximized"]			= self.isMaximized()
		self.config["layout_window_width"]		= self.geometry().width()
		self.config["layout_window_height"]		= self.geometry().height()
		self.config["layout_xposition"]			= self.pos().x()
		self.config["layout_yposition"]			= self.pos().y()

		self.config["files_reload"]				= self.save_files_action.isChecked()
		self.config["files_experimental"]		= json.dumps(self.config["files_experimental"])
		self.config["files_prediction"]			= json.dumps(self.config["files_prediction"])
		self.config["files_assigned"]			= json.dumps(self.config["files_assigned"])
		self.config["files_table"]				= self.main_widget._catalogue.to_json(orient="values")

		output_dict = {}
		for key, value in self.config.items():
			category, name = key.split("_", 1)
			category = category.capitalize()
			if category not in output_dict:
				output_dict[category] = {}
			output_dict[category][name] = value

		if self.config["files_reload"] != True:
			del output_dict["Files"]

		config = configparser.ConfigParser()
		for section in output_dict:
			config.add_section(section)
			for key in output_dict[section]:
				config.set(section, key, str(output_dict[section][key]))

		with open(f"{app_tag}.ini", "w+", encoding = "utf-8") as file:
			config.write(file)
			
		self.main_widget.notification("Options were saved successfully!")

	def load_options(self):
		self.config = {}

		def get_config(config, name, default=None, numeric=False, integer=False, json_flag=False, color=False, colors=False):
			category, name = name.split("_", 1)
			category = category.capitalize()
			if category in config:
				if name in config[category]:
					val = config[category][name]
					if val == "True":
						val = True
					elif val == "False":
						val = False
					elif val == "None":
						val = None
					if numeric == True or integer == True:
						try:
							if integer == True:
								val = int(val)
							else:
								val = np.float64(val)
						except ValueError:
							if default != None:
								val = default
							else:
								val = 0
					elif json_flag == True:
						try:
							val = json.loads(val)
						except json.decoder.JSONDecodeError:
							val = default
						except TypeError:
							val = default
					elif color == True:
						if not re.search(r'^#(?:[0-9a-fA-F]{3}){1,2}$', val):
							val = default
					elif colors == True:
						try:
							tmp_colors = json.loads(val)
							tmp_colors = [color for color in tmp_colors if re.search(r'^#(?:[0-9a-fA-F]{3}){1,2}$', color)]
							if len(tmp_colors) > 0:
								val = tmp_colors
							else:
								val = default
						except:
							val = default
					return(val)
			return(default)

		try:
			self.__suspend_updating__ = True
			config = configparser.ConfigParser()
			config.read(f"{app_tag}.ini")

			config_specs = {
				"layout_direction" :				["horizontal"],
				"layout_options_boxes" :			[True],
				"layout_lines_box" :				[True],
				"layout_log_box" :					[True],
				"layout_style" :					["light"],
				"layout_window_width" :				[None, "integer"],
				"layout_window_height" :			[None, "integer"],
				"layout_maximized" :				[False],
				"layout_xposition" :				[None, "integer"],
				"layout_yposition" :				[None, "integer"],
				"layout_scheme" :					[None],
				
				"color_experimental" :				["#000000", "color"],
				"color_assigned" :					["#ff38fc", "color"],
				"color_current" :					["#71eb34", "color"],
				"color_catfiles" :					[["#ff1a29", "#e6254b", "#de3aa2", "#de763a", "#ded03a"], "colors"],
				"color_fitcolor" :					["#fab520", "color"],
				
				"fit_error" :						[1, "integer"],
				"fit_engine" :						["polynom"],
				"fit_alwaysfit" :					[False],
				"fit_errorfactor" :					[0.001, "numeric"],
				"fit_rankpolynom" :					[5, "integer"],
				"fit_alwaysassign" :					[True],
				
				"plot_dpi" :						[100, "numeric"],
				"plot_annotations_dict" :			[{"x": 1, "y": 1, "horizontalalignment": "right", "verticalalignment": "top"}, "json_flag"],
				"plot_font_dict" :					[{"size":10}, "json_flag"],
				"plot_matplotlibkwargs" :			[{},"json_flag"],
				"plot_mark_assigned" :				[True],
				"plot_yrange_exp" :					[2, "integer"],
				"plot_yrange_lin" :					[2, "integer"],
				"plot_show_xaxis" :					[False],
				"plot_coupled" :					[True],
				"plot_show_yaxis_exp" :				[False],
				"plot_show_yaxis_lin" :				[False],
				"plot_ymargin" :					[0.1, "numeric"],
				"plot_annotate" :					[True],
				"plot_hover" :						[True],
				"plot_hover_cutoff" :				[20, "numeric"],
				"plot_number_of_plots" :			[1, "integer"],
				"plot_tight_layout" :				[False],
				"plot_expyzoom" :					[0, "numeric"],
				"plot_linyzoom" :					[0, "numeric"],
				"plot_hidecatalogue" :				[False],

				"series_number_of_plots" :			[1, "integer"],
				"series_width" :					[100, "integer"],
				"series_ignore" :					[False],
				"series_qns" :						[3, "integer"],
				"series_J0" :						[0, "integer"],
				"series_expression" :				[None],
				"series_transition" :				[[[1, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0], [True, False, True, False, False, False]], "json_flag"],
				"series_method" :					[1, "integer"],

				"flag_append" :						[True],
				"flag_automatic_draw" :				[True],
				"flag_binning" :					[4000, "integer"],
				"flag_xcolumn" :					[0, "integer"],
				"flag_ycolumn" :					[1, "integer"],
				"flag_separator" :					[9, "integer"],
				"flag_reverse_shift" :				[False],
				"flag_time_plots" :					[False],
				"flag_debug" :						[False],
				"flag_alwaysshowlog" :				[True],

				"pipe_command" :					[None],
				"pipe_reread" :						[True],
				
				"blendedlineswindow_lineshape" :	["Gauss"],
				"blendedlineswindow_derivative" :	[0, "integer"],
				"blendedlineswindow_transparency" :	[0.2, "numeric"],
				"blendedlineswindow_maxfwhm" :		[""],
				"blendedlineswindow_polynom" :		[0, "integer"],
				"blendedlineswindow_showbaseline" :	[True],
				
				"lineshapewindow_lineshape" :		["Gauss"],
				"lineshapewindow_derivative" :		[0, "integer"],
				"lineshapewindow_gauss" :			[""],
				"lineshapewindow_lorentz" :			[""],
				
				"pathfinderwindow_start" :			[""],
				"pathfinderwindow_stop" :			[""],
				"pathfinderwindow_results" :		[10, "integer"],
				"pathfinderwindow_condition" :		[""],
				"pathfinderwindow_types" :			[json.dumps(["false"]*3)],
				
				"seriesfitwindow_qns" :				[3, "integer"],
				"seriesfitwindow_transition" :		[None, "json_flag"],
				"seriesfitwindow_function" :		[""],
				
				"files_reload" :					[True],
				"files_experimental" :				[{}, "json_flag"],
				"files_prediction" :				[{}, "json_flag"],
				"files_assigned" :					[{}, "json_flag"],
				"files_table" :						[None],
			}

			for key, params in config_specs.items():
				args_dict = {"default" : params[0]}
				if len(params) > 1:
					args_dict[params[1]] = True
				self.config[key] = get_config(config, key, **args_dict)

			## Initialize Plot Widget
			##
			self.main_widget = PlotWidget(self)
			self.setCentralWidget(self.main_widget)
			self.create_menu()
			## Plot
			matplotlib.rc('font', **self.config["plot_font_dict"])
			self.mark_assigned_action.setChecked(self.config["plot_mark_assigned"])
			self.main_widget.radio_yexp.buttons()[self.config["plot_yrange_exp"]].setChecked(True)
			self.main_widget.radio_ylin.buttons()[self.config["plot_yrange_lin"]].setChecked(True)
			self.main_widget.checkbox_xaxis.setChecked(self.config["plot_show_xaxis"])
			self.main_widget.checkbox_coupled.setChecked(self.config["plot_coupled"])
			self.main_widget.checkbox_show_yaxis_exp.setChecked(self.config["plot_show_yaxis_exp"])
			self.main_widget.checkbox_show_yaxis_lin.setChecked(self.config["plot_show_yaxis_lin"])
			self.main_widget.checkbox_annotate.setChecked(self.config["plot_annotate"])
			self.main_widget.checkbox_hover.setChecked(self.config["plot_hover"])
			self.main_widget.alter_plots(absolute=self.config["plot_number_of_plots"], update = False)
			self.main_widget.expScaleZoomInput.setValue(self.config["plot_expyzoom"])
			self.main_widget.linScaleZoomInput.setValue(self.config["plot_linyzoom"])
			## Layout
			self.main_widget.topRightGroupBox.setHidden(not self.config["layout_options_boxes"])
			self.main_widget.topLeftGroupBox.setHidden(not self.config["layout_options_boxes"])
			self.main_widget.catalogueBox.setHidden(not self.config["layout_lines_box"])
			self.main_widget.catalogueBox.setHidden(not self.config["layout_log_box"])
			self.change_style(self.config["layout_style"])
			## Series
			self.main_widget.alter_QNs(0)
			self.main_widget.series_nop_input.setValue(self.config["series_number_of_plots"])
			self.main_widget.series_width_input.setValue(self.config["series_width"])
			self.main_widget.series_ignore_checkbox.setChecked(self.config["series_ignore"])
			self.main_widget.input_J0.setValue(self.config["series_J0"])
			self.main_widget.inputExpr.setText(self.config["series_expression"])
			self.main_widget.tabs.setCurrentIndex(self.config["series_method"])
			if self.config["series_transition"] != None:
				[inp.setValue(val) for val, inp in zip(self.config["series_transition"][0], self.main_widget.QNUs)]
				[inp.setValue(val) for val, inp in zip(self.config["series_transition"][1], self.main_widget.QNLs)]
				[cb.setChecked(val) for val, cb in zip(self.config["series_transition"][2], self.main_widget.QNcbs)]
			## Fit
			self.main_widget.error_field.setValue(self.config["fit_error"])
			self.always_fit_action.setChecked(self.config["fit_alwaysfit"])
			self.always_assign_action.setChecked(self.config["fit_alwaysassign"])
			## Flags
			self.save_lines_append_action.setChecked(self.config["flag_append"])
			self.automatic_draw_action.setChecked(self.config["flag_automatic_draw"])
			self.main_widget.binning_field.setValue(self.config["flag_binning"])
			self.always_show_log_action.setChecked(self.config["flag_alwaysshowlog"])
			## Files
			if self.config["files_reload"] == True:
				self.save_files_action.setChecked(True)
				should_reread = False
				if self.config["files_experimental"] != None:
					should_reread = True
					
				if self.config["files_prediction"] != None:
					should_reread = True

				if self.config["files_assigned"] != None:
					should_reread = True
				
				if self.config["files_table"] != None:
					tmp_df = pd.read_json(self.config["files_table"], orient="values")
					for i, row in tmp_df.iterrows():
						self.main_widget._catalogue.loc[i] = row.values
					self.main_widget.catalogueModel.update()
				if should_reread == True:
					self.reread_files(do_QNs=True)
			## Layout of Window
			if self.config["layout_maximized"] == True:
				self.showMaximized()
			elif self.config["layout_window_width"] != None and self.config["layout_window_height"] != None:
				self.resize(self.config["layout_window_width"], self.config["layout_window_height"])

			if self.config["layout_xposition"] != None and self.config["layout_yposition"] != None:
				self.move(self.config["layout_xposition"], self.config["layout_yposition"])
			else:
				self.move(0,0)
		except Exception as E:
			self.main_widget.notification("Your option file seems to be corrupted, as there occurred an error when loading the specified settings. Please repair or delete the option.ini file or save your defaults again now.")
			exc_type, exc_value, exc_traceback = sys.exc_info()
			self.main_widget.notification(f"The exception was '{str(E)}' and the traceback was: \n{''.join(tb.format_tb(exc_traceback))}")
			
		finally:
			self.__suspend_updating__ = False

	@threading_d
	@working_d
	@synchronized(locks["exp_df"])
	def load_exp(self, keep_old = False, add_files = False, reread = False, skip_update=False):
		if reread == True:
			keep_old = False
			fnames = self.config["files_experimental"].keys()
		elif add_files != False:
			fnames = add_files
		else:
			fnames = []
			
		if keep_old == False:
			self.main_widget._exp_data_df.drop(self.main_widget._exp_data_df.index, inplace = True)
			if reread == False:
				self.config["files_experimental"].clear()
		
		sep = chr(self.config["flag_separator"])
		for fname in fnames:
			try:
				df = self.main_widget._exp_data_df.copy()
				files = self.config["files_experimental"].copy()
				
				if fname in files:
					df.drop(df[df.filename==fname].index, inplace=True)
				else:
					files[fname] = {}
				
				data = pd.read_csv(fname, sep=sep, dtype=np.float64, engine="c")
				if len(data.columns) < 2:
					self.main_widget.notification(f"The file {fname} has less than two columns, it was skipped.")
					continue
				column_names = [i for i in range(len(data.columns))]
				if self.config["flag_xcolumn"] in range(len(data.columns)):
					column_names[self.config["flag_xcolumn"]] = "x"
				else:
					self.main_widget.notification(f"The specified value for the x_column is higher than the number of columns in {fname}, therefore the first column will be used.")
					column_names[0] = "x"
				if self.config["flag_ycolumn"] in range(len(data.columns)):
					column_names[self.config["flag_ycolumn"]] = "y"
				else:
					self.main_widget.notification(f"The specified value for the y_column is higher than the number of columns in {fname}, therefore the second column will be used.")
					column_names[1] = "y"
				data.columns = column_names
				if files[fname].get("invert", False) == True:
					data["y"] = -data["y"]
				data["filename"] = fname
				df = df.append(data[["x","y","filename"]], ignore_index = True)
				
				if files[fname].get("color") == None:
					files[fname]["color"] = self.config["color_experimental"]
			except Exception as E:
				df = self.main_widget._exp_data_df.copy()
				files = self.config["files_experimental"].copy()
				self.main_widget.notification(f"There occurred an error when loading the spectrum from {fname}. Please check the file.")
				if self.config["flag_debug"] == True:
					self.main_widget.notification(f"The error message reads: {str(E)}")
				continue

			self.main_widget._exp_data_df = df
			self.config["files_experimental"].clear()
			self.config["files_experimental"].update(files)
		
		df = self.main_widget._exp_data_df
		df.drop_duplicates(keep = "first", inplace = True)
		df.sort_values("x", inplace = True)

		#Preprocess maxima
		y_min = df["y"].min()
		y_max = df["y"].max()
		self.main_widget._yrange = [y_min, y_max]
		self.main_widget.signalclass.updateWindows.emit()
		if skip_update != True:
			self.main_widget.update_plots()

	@threading_d
	@working_d
	@synchronized(locks["lin_df"])
	def load_cat(self, keep_old = False, add_files = False, do_QNs = True, reread = False, skip_update=False):
		def letter_to_number(val):
			if val == "":
				val = -1
			elif val[0].isalpha():
				val = str(ord(val[0].upper())-55)+val[1:]
			return(np.int16(val))

		if reread == True:
			keep_old = False
			fnames = self.config["files_prediction"].keys()
		elif add_files != False:
			fnames = add_files
		else:
			fnames = []
		
		if keep_old == False:
			self.main_widget._lin_data_df.drop(self.main_widget._lin_data_df.index, inplace=True)
			if reread == False:
				self.config["files_prediction"].clear()
		

		dtypes = self.main_widget._lin_data_df_dtypes
		widths = [13, 8, 8, 2, 10, 3, 7, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
		column_names = self.main_widget._lin_data_df_columns
		preset_colors = self.config["color_catfiles"]
		dtypes_dict = {column_names[i]:dtypes[i] for i in range(len(column_names))}
		converters = {key:letter_to_number for key in column_names[8:20]}
		for key in column_names[8:20]:
			del dtypes_dict[key]
		for fname in fnames:
			try:
				df = self.main_widget._lin_data_df.copy()
				files = self.config["files_prediction"].copy()
				
				if fname in files:
					df.drop(df[df.filename==fname].index, inplace=True)
				else:
					files[fname] = {}
				
				data = pd.read_fwf(fname, widths=widths, names=column_names, converters=converters, skip_blank_lines=True, dtype=dtypes_dict)
				data["filename"] = fname
				data["y"] = 10 ** data["y"]
				df = df.append(data, ignore_index = True)
				
				if files[fname].get("color") == None:
					files[fname]["color"] =  preset_colors[list(files).index(fname)%len(preset_colors)]
			except Exception as E:
				df = self.main_widget._lin_data_df.copy()
				files = self.config["files_prediction"].copy()
				self.main_widget.notification(f"There occurred an error when loading the lines from {fname}. Please check the file.")
				if self.config["flag_debug"] == True:
					self.main_widget.notification(f"The error message reads: {str(E)}")
				continue
			
			self.main_widget._lin_data_df = df
			self.config["files_prediction"].clear
			self.config["files_prediction"].update(files)
			
		df = self.main_widget._lin_data_df
		df.drop_duplicates(keep = "first", inplace = True)
		df.sort_values("x", inplace = True)
		
		if do_QNs == True and len(df) != 0:
			qn_labels = ['qnu1', 'qnu2', 'qnu3', 'qnu4', 'qnu5', 'qnu6']
			QNs = len(qn_labels)
			for i in range(len(qn_labels)):
				tmp_label = qn_labels[i]
				tmp_unique = df[tmp_label].unique()
				if len(tmp_unique) == 1 and tmp_unique[0] == -1:
					QNs = i
					break
			self.main_widget.alter_QNs(absolute=QNs)
			self.main_widget.notification(f"After analysing your cat files the number of QNs was set to {QNs}.")

		# Preprocessing Maxima and No binning case
		y_min = df["y"].min()
		y_max = df["y"].max()
		self.main_widget._tyrange = [y_min, y_max]
		self.main_widget.signalclass.updateWindows.emit()
		if skip_update != True:
			self.main_widget.update_plots()

	@threading_d
	@working_d
	@synchronized(locks["ass_df"])
	def load_lin(self, keep_old = False, add_files = False, reread = False, skip_update=False):
		def column_to_float(val):
			val = val.strip()
			if val == "":
				return(-1)
			elif val[0].isalpha():
				val = str(ord(val[0].upper())-55)+val[1:]
			return(np.float64(val))

		if reread == True:
			keep_old = False
			fnames = self.config["files_assigned"].keys()
		elif add_files != False:
			fnames = add_files
		else:
			fnames = []
			
		if keep_old == False:
			self.main_widget._ass_data_df.drop(self.main_widget._ass_data_df.index, inplace = True)
			if reread == False:
				self.config["files_assigned"].clear()
		
		dtypes = self.main_widget._ass_data_df_dtypes
		catindices = range(0,37,3)
		column_names = self.main_widget._ass_data_df_columns
		dtypes_dict = {column_names[i]:dtypes[i] for i in range(len(column_names))}
		for fname in fnames:
			try:
				df = self.main_widget._ass_data_df.copy()
				files = self.config["files_assigned"].copy()
				
				if fname in files:
					df.drop(df[df.filename==fname].index, inplace=True)
				else:
					files[fname] = {}
					
				data = []
				with open(fname, "r") as file:
					for line in file:
						if line.strip() == "" or line.startswith("#"):
							continue
						tmp_line_content = [column_to_float(line[i:j]) for i,j in zip(catindices[:-1], catindices[1:])]+[column_to_float(x) for x in line[36:].split()]
						tmp_line_content.extend([-1]*(15-len(tmp_line_content)))
						data.append(tmp_line_content)
				data = pd.DataFrame(data, dtype = np.float64)
				data["filename"] = fname
				data.columns = column_names
				data = data.astype(dtypes_dict)
				## Set correct columns for data
				qn_labels = column_names[0:12]
				QNs = len(qn_labels)
				for i in range(len(qn_labels)):
					tmp_unique = data[qn_labels[i]].unique()
					if len(tmp_unique) == 1 and tmp_unique[0] == -1:
						QNs = i
						break
				QNs = int(QNs/2)
				columns_qn = [f"qnu{i+1}" for i in range(QNs)]+[f"qnl{i+1}" for i in range(QNs)]+[f"qnu{i+1}" for i in range(QNs, 6)]+[f"qnl{i+1}" for i in range(QNs, 6)]
				new_column_names = columns_qn + column_names[12:]
				data.columns = new_column_names
				df = df.append(data, ignore_index = True, sort=False)
				
				if files[fname].get("color") == None:
					files[fname]["color"] = self.config["color_assigned"]
			except Exception as E:
				df = self.main_widget._ass_data_df.copy()
				files = self.config["files_assigned"].copy()
				self.main_widget.notification(f"There occurred an error when loading the lines from {fname}. Please check the file.")
				if self.config["flag_debug"] == True:
					self.main_widget.notification(f"The error message reads: {str(E)}")
				continue
				
			self.main_widget._ass_data_df = df	
			self.config["files_assigned"].clear()
			self.config["files_assigned"].update(files)
	
		df = self.main_widget._ass_data_df
		df.drop_duplicates(keep = "first", inplace = True)
		df.sort_values("x", inplace = True)

		self.main_widget.signalclass.updateWindows.emit()
		if skip_update != True:
			self.main_widget.update_plots()

	def reread_files(self, do_QNs=False):
		kwargs = {"reread": True, "skip_update":True}
		threads = []
		threads.append(self.load_exp(**kwargs))
		threads.append(self.load_cat(do_QNs=do_QNs, **kwargs))
		threads.append(self.load_lin(**kwargs))
		if self.__suspend_updating__ == True:
			for thread in threads:
				thread.join()
		self.main_widget.update_plots()
	
	def change_style(self, style = None):
		styles = ["light", "dark", "custom"]
		if style == None:
			self.config["layout_style"] = styles[(styles.index(self.config["layout_style"])+1)%len(styles)]
		elif style in styles:
			self.config["layout_style"] = style
		else:
			self.config["layout_style"] = styles[0]
		if self.config["layout_scheme"] == None and self.config["layout_style"] == "custom":
			self.config["layout_style"] = "light"
			
		if self.config["layout_style"] == "light":
			app.setStyle("WindowsVista")
			app.setPalette(app.style().standardPalette())
			matplotlib.style.use("default")
			self.main_widget.fig.patch.set_facecolor('white')
			self.main_widget.alter_plots(absolute = len(self.main_widget._taxs))

		elif self.config["layout_style"] in ["custom", "dark"]:
			if self.config["layout_style"] == "dark":
				colors = {
					"color_bg": 		QColor(53, 53, 53),
					"color_fg":			Qt.white,
					"color_con_bg":		QColor(25, 25, 25),
					"color_bright":		Qt.red,
					"color_highlight":	QColor("#5a65ec"),
					"color_con_active":	Qt.black,
					"matplotlib_fc":	"black",
					"matplotlib_style":	"dark_background",
				}
			else:
				colors = {k: QColor(v) if k.startswith("color") else v for k, v in json.loads(self.config["layout_scheme"].items())}
			dark_palette = QPalette()
			dark_palette.setColor(QPalette.Window, colors["color_bg"])
			dark_palette.setColor(QPalette.WindowText, colors["color_fg"])
			dark_palette.setColor(QPalette.Base, colors["color_con_bg"])
			dark_palette.setColor(QPalette.AlternateBase, colors["color_bg"])
			dark_palette.setColor(QPalette.ToolTipBase, colors["color_fg"])
			dark_palette.setColor(QPalette.ToolTipText, colors["color_bg"])
			dark_palette.setColor(QPalette.Text, colors["color_fg"])
			dark_palette.setColor(QPalette.Button, colors["color_bg"])
			dark_palette.setColor(QPalette.ButtonText, colors["color_fg"])
			dark_palette.setColor(QPalette.BrightText, colors["color_bright"])
			dark_palette.setColor(QPalette.Link, colors["color_highlight"])
			dark_palette.setColor(QPalette.Highlight, colors["color_highlight"])
			dark_palette.setColor(QPalette.HighlightedText, colors["color_con_active"])

			app.setStyle("Fusion")
			app.setPalette(dark_palette)
			matplotlib.style.use(colors["matplotlib_style"])
			self.main_widget.fig.patch.set_facecolor(colors["matplotlib_fc"])
			self.main_widget.alter_plots(absolute = len(self.main_widget._taxs))

	def go_to(self):
		resp, rc = QInputDialog.getText(self, 'Go to frequency', 'Frequency:')
		if rc == True:
			if self.main_widget._lastclickedplot not in self.main_widget._axs.keys():
				 self.main_widget._lastclickedplot = 0
			if self.main_widget.checkbox_coupled.isChecked() == False:
				offset = self.main_widget.offset(self.main_widget._lastclickedplot)
			else:
				offset = self.main_widget.offset(len(self.main_widget._axs)-1)
			try:
				resp = float(resp)
				resp = resp - offset
				self.main_widget.update_plots({"step":{"abs":resp}})
			except ValueError:
				pass

	def set_zoom(self):
		resp, rc = QInputDialog.getText(self, 'Set view width', 'Width:')
		if rc == True:
			try:
				resp = float(resp)
				self.main_widget.update_plots({"zoom":{"abs":resp}})
			except ValueError:
				pass

	def plot_number(self):
		resp, rc = QInputDialog.getText(self, 'How many plots do you want: ', 'Number:')
		if rc == True:
			try:
				resp = float(resp)
				self.main_widget.alter_plots(absolute=resp)
			except ValueError:
				pass

	@synchronized(locks["cat_df"])
	def save_lines_cat(self, path=None, force_new=False):
		append = self.save_lines_append_action.isChecked()
		if append:
			options = {"options" : QFileDialog.DontConfirmOverwrite}
		else:
			options = {}
		if path == None:
			path, ext = QFileDialog.getSaveFileName(self, 'Save file', '', **options)
			if path == "":
				return
			path = path.encode('utf-8')

		catalogue = self.main_widget._catalogue
		output = ""

		if append == True and force_new != True:
			handle = "a+"
		else:
			handle = "w+"

		for index, row in catalogue.iterrows():
			freq = row["x"]
			if freq > 99999999.9999:
				freq = 99999999.9999
			elif freq < 0.0001:
				freq = 0.0001
			error = row["error"]
			if error > 999.9999:
				error = 999.9999
			elif error < -99.9999:
				error = -99.9999
			intens = np.log10(row["y"]) if row["y"] > 0 else 0
			if intens < -99.999:
				intens = -99.999
			elif intens > 999.999:
				intens = 999.999

			QNs_string = ""
			for QNlabel in ['qnu1', 'qnu2', 'qnu3', 'qnu4', 'qnu5', 'qnu6', 'qnl1', 'qnl2', 'qnl3', 'qnl4', 'qnl5', 'qnl6']:
				if row[QNlabel] in ["", -1]:
					QNs_string += "  "
				else:
					QNs_string += f"{row[QNlabel]:2.0f}"

			output += f"{freq:13.4f}{error:8.4f}{intens:8.4f}{row['degfreed']:2.0f}{row['elower']:10.4f}{row['usd']:3.0f}{row['tag']:7.0f}{row['qnfmt']:4.0f}{QNs_string}"
			output += "\n"
		with open(path, handle) as file:
			file.write(output)

	@synchronized(locks["cat_df"])
	def save_lines_lin(self, path = None):
		append = self.save_lines_append_action.isChecked()
		if append:
			options = {"options" : QFileDialog.DontConfirmOverwrite}
		else:
			options = {}
		if path == None:
			path, ext = QFileDialog.getSaveFileName(self, 'Save file', '', **options)
			if path == "":
				return
			path = path.encode('utf-8')
		catalogue = self.main_widget._catalogue
		output = ""

		if append == True:
			handle = "a+"
		else:
			handle = "w+"

		for index, row in catalogue.iterrows():
			QNs_string = ""
			pad_string = ""
			for QNlabel in ['qnu1', 'qnu2', 'qnu3', 'qnu4', 'qnu5', 'qnu6', 'qnl1', 'qnl2', 'qnl3', 'qnl4', 'qnl5', 'qnl6']:
				if row[QNlabel] in ["", -1]:
					pad_string += "   "
				else:
					QNs_string += f"{row[QNlabel]:3.0f}"
			QNs_string = QNs_string + pad_string
			output += f"{QNs_string} {row['x']:13.4f} {row['error']:8.4f}"
			output += "\n"
		with open(path, handle) as file:
			file.write(output)

	@synchronized(locks["cat_df"])
	def clear_lines(self):
		self.main_widget._catalogue.drop(self.main_widget._catalogue.index, inplace = True)
		self.main_widget.catalogueModel.update()
		self.main_widget.update_plots()

	def change_fitengine(self, new = None):
		available = self.main_widget.fitfunctions
		current = self.config["fit_engine"]
		if current not in available:
			current = available[0]
		if new != None and new in available:
			new_index = available.index(new)
		else:
			new_index = (available.index(current)+1)%len(available)
		self.config["fit_engine"] = available[new_index]
		self.main_widget.notification(f"Fitting with {available[new_index]}")
		
	def change_fitcolor(self):
		color = QColorDialog.getColor()
		color = match_color(color.name())
		if color != False:
			self.config["color_fitcolor"] = color			

	@synchronized(locks["windows"])
	def open_window(self, windowname):
		if self.open_windows.get(windowname) == None:
			self.open_windows[windowname] = self.window_class.get(windowname)(windowname)
		self.open_windows[windowname].show()
		self.open_windows[windowname].activateWindow()

	@synchronized(locks["pipe"])
	def run_pipe(self):
		command = self.config["pipe_command"]
		if type(command) == str:
			command = command.replace("\n", " && ")
			try:
				output = subprocess.check_output(command, shell = True)
				output = output.decode("utf-8")

				self.main_widget.notification(f"The subprocess was started and returned:\n{output}")
				if self.config.get("pipe_reread") == True:
					self.reread_files()
			except Exception as E:
				message.append(f"The subprocess failed with the Exception '{str(E)}'.")
				message.append("Command '{}' returned with error (code {}) and output: {}".format(E.cmd, E.returncode, E.output))
				self.main_widget.notification("\n".join(message))
		else:
			self.main_widget.notification("No Pipe command specified, as a result no Pipe process was started.")

	def send_mail_to_author(self):
		quotes = json.loads(quotes_str)
		quote = quotes[random.randint(0,len(quotes)-1)]
		quote = quote.replace("\n", "%0d%0a")
		webbrowser.open(f"mailto:bonah@ph1.uni-koeln.de?subject=LWP: &body=%0d%0a%0d%0a%0d%0a%0d%0a{quote}")


	##
	## Menu Setup and Functions Functions
	##
	def create_menu(self):
		self.file_menu = self.menuBar().addMenu("&Files")
		self.view_menu = self.menuBar().addMenu("&View")
		self.line_menu = self.menuBar().addMenu("&Lines")
		self.fit_menu = self.menuBar().addMenu("&Fit")
		self.plot_menu = self.menuBar().addMenu("&Plot")
		self.modules_menu = self.menuBar().addMenu("&Modules")
		self.info_menu = self.menuBar().addMenu("&Info")

		##
		## File Menu
		load_exp_action = self.create_action("&Load Spectrum", shortcut="Ctrl+O", slot=self.load_exp, tip="Open a spectrum")
		add_file_action = self.create_action("&Add Spectrum", shortcut="Shift+O", slot=lambda x:self.load_exp(keep_old = True), tip="Add a spectrum")
		load_cate_action = self.create_action("&Load Cat File", shortcut="Ctrl+L", slot=self.load_cat, tip="Open a cat file")
		add_line_action = self.create_action("&Add Cat File", shortcut="Shift+L",slot=lambda x:self.load_cat(keep_old=True), tip="Add a cat file")
		load_lin_action = self.create_action("&Load Lin File", shortcut="Ctrl+I", slot=self.load_lin, tip="Load lin file")
		add_assigned_action = self.create_action("&Add Lin File", shortcut="Shift+I", slot=lambda x:self.load_lin(keep_old=True), tip="Load already assigned lines")
		reread_files_action = self.create_action("&Reread Files", slot = self.reread_files, tip="Reread all external files", shortcut="Shift+R")
		edit_exp_action = self.create_action("&Edit Spectrum Files", shortcut="Shift+5", tip="Press to see list of used spectrum files", slot=lambda a, x="spectraDialog_window": self.open_window(x))
		edit_lines_action = self.create_action("&Edit Cat Files", shortcut="Shift+6", tip="Change colors of or delete cat files", slot=lambda a, x="linesDialog_window": self.open_window(x))
		edit_ass_action = self.create_action("&Edit Lin Files", shortcut="Shift+7", tip="Press to see list of used lin files", slot=lambda a, x="assignedDialog_window": self.open_window(x))
		self.save_files_action = self.create_action("&Save Files", tip="Save current files as options", checkable=True)
		save_options_action = self.create_action("&Save Options", shortcut="Ctrl+D", tip="Save view options as default", slot=lambda a:self.save_options())
		quit_action = self.create_action("&Quit", slot=self.close, tip="Close the application")

		##
		## View Menu
		toggle_toolbar_action = self.create_action("&Toggle MPL Toolbar", shortcut="Shift+1", tip="Show or hide toolbar to edit or save plot", slot=lambda a, x=self.main_widget.mpl_toolbar: x.setHidden(not x.isHidden()))
		toggle_options_action = self.create_action("&Toggle Options Boxes", shortcut="Shift+2", tip="Toggle visibility for options boxes.", slot=lambda a,x=self.main_widget.topLeftGroupBox, y=self.main_widget.topRightGroupBox: [elem.setHidden(not elem.isHidden()) for elem in (x,y)])
		toggle_lines_action = self.create_action("&Toggle Lines Table", shortcut="Shift+3", tip="Toggle visibility of assigned lines table", slot=lambda a,x=self.main_widget.catalogueBox: x.setHidden(not x.isHidden()))
		toggle_log_action = self.create_action("&Toggle Log", shortcut="Shift+4", tip="Toggle visibility of log tab", slot=lambda a,x=self.main_widget.logBox: x.setHidden(not x.isHidden()))
		change_style_action = self.create_action("&Change Style", tip="What we know is a drop. What we do not know… is an ocean.", slot=lambda a:self.change_style())
		self.always_show_log_action = self.create_action("&Force Show Log", tip="Show log window on new message", checkable= True)

		##
		## Lines Menu
		save_lines_cat_action = self.create_action("&Save Lines Cat", shortcut="Ctrl+E", tip="Save the assigned lines to a cat file.", slot=lambda a:self.save_lines_cat())
		save_lines_lin_action = self.create_action("&Save Lines Lin", shortcut="Ctrl+R", tip="Save the assigned lines to a lin file.", slot=lambda a:self.save_lines_lin())
		self.save_lines_append_action = self.create_action("&Append", tip="Append to file instead of overwriting.", checkable = True)
		clear_lines_action = self.create_action("&Clear Lines", tip="Clear the assigned lines.", slot=self.clear_lines)


		##
		## Fit Menu
		change_fitengine_action = self.create_action("&Change Function", shortcut="Shift+F", tip="Change the shape of the underlying function", slot=lambda x:self.change_fitengine())
		self.always_fit_action = self.create_action("&Always Fit", tip="Always fit on selecting a range, helps if you have troubles with modifiers", checkable = True)
		change_fit_color_action = self.create_action("&Change Fit Color", slot=lambda x:self.change_fitcolor())
		self.always_assign_action = self.create_action("&Always Assign", tip="Always assign fitted line", checkable = True)
		manual_assign_action = self.create_action("&Manual Assign", tip="Manual assign when always assign is off", slot=lambda x:self.main_widget.manual_assign(), shortcut="Ctrl+Return")

		##
		## Plot Menu
		plot_number_action = self.create_action("&# Plots", shortcut="Ctrl+N", tip="Change amount of plots",slot=self.plot_number)
		go_to_action = self.create_action("&Go to ..", tip="Go to specific center frequency", shortcut="Ctrl+G",slot=self.go_to)
		set_zoom_action = self.create_action("&Set Zoom", tip="Set specific width", shortcut="Ctrl+W",slot=self.set_zoom)
		self.automatic_draw_action = self.create_action("&Automatic Draw", tip="Draw canvas automatically.", checkable = True)
		draw_manually_action = self.create_action("&Manual Draw", tip="Draw canvas.", slot=lambda a: self.main_widget.update_plots({"Force_Draw":True}))
		self.mark_assigned_action = self.create_action("&Mark Assigned", tip="Show marker when line is assigned.", checkable = True)
		lineshape_action = self.create_action("&Show Lineshape", shortcut="Shift+8", tip="Press to see lineshape", slot=lambda a, x="lineshape_window": self.open_window(x))
		blendedlines_action = self.create_action("&Blended Lines", shortcut="Shift+B", tip="Press to check for blend", slot=lambda a, x="blendedlines_window": self.open_window(x))

		##
		## Pathfinder Menu
		pathfinder_action = self.create_action("&Show Pathfinder", shortcut="Shift+9", tip="Open pathfinder window", slot=lambda a, x="pathfinder_window": self.open_window(x))

		##
		## Pipes Menu
		option_pipe_action = self.create_action("&Set up Pipe", tip="Define which action to do", slot=lambda a, x="pipe_window": self.open_window(x))
		run_pipe_action = self.create_action("&Run Pipe", tip="Start Pipe action", slot=self.run_pipe, shortcut = "Ctrl+P")
		
		##
		## Series Fit Menu
		series_fit_action = self.create_action("&Series Fit", tip="Open Series Fit Window", slot=lambda a:self.open_window("series_fit_window"))


		##
		## Contact Menu
		self.send_mail_action = self.create_action("&Send Mail to Author", tip="Send a Mail for all the good and bad things", slot=lambda x: self.send_mail_to_author())
		self.credit_action = self.create_action("&Credits", tip="See how it was made possible", slot=lambda a, x="credits_window": self.open_window(x))


		self.add_actions(self.file_menu, (load_exp_action, add_file_action, load_cate_action, add_line_action, None, load_lin_action, add_assigned_action, None, reread_files_action, None, edit_exp_action,edit_lines_action, edit_ass_action, None, self.save_files_action, save_options_action, None, quit_action))
		self.add_actions(self.view_menu, (toggle_toolbar_action, toggle_options_action, toggle_lines_action, toggle_log_action, None, change_style_action, None, self.always_show_log_action))
		self.add_actions(self.line_menu, (save_lines_cat_action, save_lines_lin_action, None, self.save_lines_append_action, None, clear_lines_action))
		self.add_actions(self.fit_menu, (self.always_fit_action, change_fitengine_action, change_fit_color_action, None, self.always_assign_action, manual_assign_action, None))
		self.add_actions(self.plot_menu, (plot_number_action, None, go_to_action, set_zoom_action, None, self.automatic_draw_action, draw_manually_action, None, self.mark_assigned_action, ))
		self.add_actions(self.modules_menu, (lineshape_action, blendedlines_action, pathfinder_action, None, option_pipe_action, run_pipe_action, None, series_fit_action))
		self.add_actions(self.info_menu, (self.send_mail_action, self.credit_action))

		choose_function_menu = QMenu("Choose Fit Function", self)
		for fn in self.main_widget.fitfunctions:
			tmp_action = self.create_action(f"&{fn}", slot=lambda a, fn=fn: self.change_fitengine(fn))
			choose_function_menu.addAction(tmp_action)

		self.fit_menu.addMenu(choose_function_menu)

	def create_action(self, text, slot=None, shortcut=None, icon=None, tip=None, checkable=False):
		action = QAction(text, self)
		if icon is not None:
			action.setIcon(QIcon(":/%s.png" % icon))
		if shortcut is not None:
			action.setShortcut(shortcut)
		if tip is not None:
			action.setToolTip(tip)
			action.setStatusTip(tip)
		if slot is not None:
			action.triggered.connect(slot)
		if checkable:
			action.setCheckable(True)
		return action

	def add_actions(self, target, actions):
		for action in actions:
			if action is None:
				target.addSeparator()
			else:
				target.addAction(action)

	##
	## Overwrite Default Class Functions
	##
	def dragEnterEvent(self, event):
		if event.mimeData().hasUrls():
			event.accept()
		else:
			event.ignore()
			
	def dropEvent(self, event):
		files = [url.toLocalFile() for url in event.mimeData().urls()]
		correct_files = []
		for file in files:
			if os.path.isfile(file):
				correct_files.append(file)
			else:
				self.main_widget.notification(f"The file {file} could not be found.")
		
		if len(correct_files) > 0:
			items = ("Spectrum","Cat","Lin")
			item, ok = QInputDialog.getItem(self, "Choose File Type","Type:", items, 0, False)
			if ok and item:
				if item == "Spectrum":
					self.load_exp(keep_old = True, add_files = correct_files)
				elif item == "Cat":
					self.load_cat(keep_old = True, add_files = correct_files)
				elif item == "Lin":
					self.load_lin(keep_old = True, add_files = correct_files)
			
	def keyPressEvent(self, event):
		self.main_widget.keyPressEvent(event)

	def closeEvent(self, event):
		if not self.main_widget._catalogue.empty:
			self.save_lines_cat(f"{app_tag}.cat", force_new=True)
		for window in list(self.open_windows.values()):
			if window != None:
				window.close()


##
## Plot Window (Main Window minus Menu Bar)
##
class PlotWidget(QWidget):
	##
	## Working Decorator
	##
	def working_d(func):
		def wrapper(self, *args, **kwargs):
			q = self.working
			ind = self.indicator
			if q.empty():
				ind.setText("Working...")
			q.put(1)
			try:
				res = func(self, *args, **kwargs)
			except Exception as E:
				exc_type, exc_value, exc_traceback = sys.exc_info()
				self.notification(f"An error occurred: {str(E)}\n {''.join(tb.format_tb(exc_traceback))}")
				res = str(E)
			finally:
				q.get()
				q.task_done()
			if q.empty():
				ind.setText("Ready")
			return(res)
		return(wrapper)
	
	def __init__(self, parent=None):
		super(PlotWidget, self).__init__(parent)

		self.main_window = parent
		self.config = self.main_window.config
		self.signalclass = SignalClass()

		self.createPlotBox()
		self.createTopLeftGroupBox()
		self.createTopRightGroupBox()
		self.createCatalogueBox()
		self.createLogBox()
		self.createQuotesBox()

		self.working = queue.Queue()
		self.currThread = None
		
		self.fitfunctions = ["polynom", "polynom_autorank", "2ndderivativeGauss", "2ndderivativeLorentz", "2ndderivativeVoigt"]

		##
		## Button Bar on Top
		##
		self.zoom_in = QPushButton("in")
		self.zoom_in.setDefault(True)
		self.zoom_in.clicked.connect(lambda x: self.update_plots({"zoom":{"factor":1/2}}))

		self.zoom_out = QPushButton("out")
		self.zoom_out.setDefault(True)
		self.zoom_out.clicked.connect(lambda x: self.update_plots({"zoom":{"factor":2}}))

		self.step_left = QPushButton("left")
		self.step_left.setDefault(True)
		self.step_left.clicked.connect(lambda x: self.update_plots({"step":{"dir":-1}}))

		self.step_right = QPushButton("right")
		self.step_right.setDefault(True)
		self.step_right.clicked.connect(lambda x: self.update_plots({"step":{"dir":+1}}))

		self.topLabel = QLabel("")
		self.topLabel.setWordWrap(True)
		self.indicator = QLabel("Ready")

		self.addPlot = QPushButton("add plot")
		self.addPlot.setDefault(True)
		self.addPlot.clicked.connect(lambda a: self.alter_plots(+1))

		self.delPlot = QPushButton("del plot")
		self.delPlot.setDefault(True)
		self.delPlot.clicked.connect(lambda a: self.alter_plots(-1))

		##
		## Top Layout
		##
		topLayout = QHBoxLayout()
		topLayout.addWidget(self.zoom_in)
		topLayout.addWidget(self.zoom_out)
		topLayout.addWidget(self.step_left)
		topLayout.addWidget(self.step_right)
		topLayout.addWidget(self.topLabel, 1)
		topLayout.addWidget(self.indicator)
		topLayout.addWidget(self.addPlot)
		topLayout.addWidget(self.delPlot)


		layout_direction = self.config["layout_direction"]
		if layout_direction == "window":
			##
			## Window Layout
			##
			mainLayout = QVBoxLayout()
			mainLayout.addLayout(topLayout)

			mainLayout.addWidget(self.plotBox)
			mainLayout.setSizeConstraint(2)
			self.setLayout(mainLayout)

		elif layout_direction == "horizontal":
			##
			## HorizontalLayout
			##
			mainLayout = QVBoxLayout()
			mainLayout.addLayout(topLayout)

			boxesLayout = QSplitter(Qt.Vertical)
			boxesLayout.setChildrenCollapsible(False)
			boxesLayout.addWidget(self.topLeftGroupBox)
			boxesLayout.addWidget(self.topRightGroupBox)
			boxesLayout.addWidget(self.catalogueBox)
			boxesLayout.addWidget(self.logBox)
			boxesLayout.addWidget(self.quotesBox)

			splitter = QSplitter(Qt.Horizontal)
			splitter.setChildrenCollapsible(False)
			splitter.addWidget(self.plotBox)
			splitter.addWidget(boxesLayout)
			splitter.setSizes((1000,100))

			mainLayout.addWidget(splitter)
			mainLayout.setSizeConstraint(2)
			self.setLayout(mainLayout)

		else:
			##
			## VerticalLayout
			##
			mainLayout = QVBoxLayout()
			mainLayout.addLayout(topLayout)

			splitterV = QSplitter(Qt.Vertical)
			splitterV.setChildrenCollapsible(False)

			horizontalBox = QSplitter(Qt.Horizontal)
			horizontalBox.setChildrenCollapsible(False)
			horizontalBox.addWidget(self.topLeftGroupBox)
			horizontalBox.addWidget(self.topRightGroupBox)

			splitterV.addWidget(self.plotBox)
			splitterV.addWidget(horizontalBox)
			splitterV.addWidget(self.catalogueBox)
			splitterV.addWidget(self.logBox)
			splitterV.addWidget(self.quotesBox)
			splitterV.setSizes((1000,100, 100, 100, 100))

			mainLayout.addWidget(splitterV)
			mainLayout.setSizeConstraint(2)
			self.setLayout(mainLayout)

	##
	## Update Canvas
	##
	def update_plots(self, dict={}):
		if self.main_window.__suspend_updating__ == True:
			return
		new_thread = threading.Thread(target = self.update_plots_worker, args=[dict])
		with locks["currThread"]:
			new_thread.start()
			self.currThread = new_thread.ident
		return(new_thread)

	@working_d
	@synchronized(locks["axs"])
	def update_plots_worker(self, dict):
		timers = []
		id = threading.current_thread().ident
		def breakpoint(self, dict, id):
			with locks["currThread"]:
				if id != self.currThread:
					raise CustomError()
		
		try:
			## Set zoom and x_0 values and alter number of plots
			timers.append(time.time())
			if "apply_series" in dict:
				if self.series_ignore_checkbox.isChecked() != True:
					self.alter_plots(absolute = int(self.series_nop_input.text()), update = False)
					width = float(self.series_width_input.text())
					dict["zoom"]={"abs": width}
			if "zoom" in dict:
				if "factor" in dict["zoom"]:
					self._range = self._range*np.float64(dict["zoom"]["factor"])
				elif "abs" in dict["zoom"]:
					abs = np.float64(dict["zoom"]["abs"])
					if not(abs > 0):
						abs = 20
					self._range = abs
			if "step" in dict:
				if "dir" in dict["step"]:
					self._x_0 = self._x_0+np.float64(dict["step"]["dir"])*self._range/2
				elif "percent" in dict["step"]:
					self._x_0 = self._x_0+np.float64(dict["step"]["percent"])/100*self._range
				elif "abs" in dict["step"]:
					self._x_0 = np.float64(dict["step"]["abs"])

			timers.append(time.time())
			breakpoint(self, dict, id)
			
			dict.update(self._queue_dict)
			self._queue_dict = {}
			
			length_axs = len(self._taxs)
			binning = self.binning_field.value()
			## If atbu is given, no x- and y-range updates are done
			axs_to_be_updated = dict.get("atbu", range(length_axs))
			self.assigned_dataframe()
			## Axis Options
			for i in axs_to_be_updated:
				curr_ax = self._axs[i]
				if self.checkbox_xaxis.isChecked() == True or i==len(self._axs)-1:
					curr_ax.get_xaxis().set_visible(True)
				else:
					curr_ax.get_xaxis().set_visible(False)
				if self.checkbox_show_yaxis_exp.isChecked() == True:
					curr_ax.get_yaxis().set_visible(True)
				else:
					curr_ax.get_yaxis().set_visible(False)
				if self.checkbox_show_yaxis_lin.isChecked() == True:
					self._taxs[i].get_yaxis().set_visible(True)
				else:
					self._taxs[i].get_yaxis().set_visible(False)

			## Layout if tight layout is off
			mplkwargs = self.config["plot_matplotlibkwargs"]
			if self.config["plot_tight_layout"] == False:
				if self.checkbox_xaxis.isChecked() == True and self._flag_plots_spaced == False:
					self.fig.subplots_adjust(wspace=0, hspace=0.3, **mplkwargs)
					self._flag_plots_spaced = True
				elif self.checkbox_xaxis.isChecked() == False and self._flag_plots_spaced == True:
					self.fig.subplots_adjust(wspace=0, hspace=0, **mplkwargs)
					self._flag_plots_spaced = False

			timers.append(time.time())
			breakpoint(self, dict, id)

			## Scale y range and y range according to settings
			x_ranges = {}
			y_ranges_exp = {}
			y_ranges_lin = {}
			
			calc_middles = self.offsets(axs_to_be_updated) + self._x_0
			for i, calc_middle in zip(axs_to_be_updated, calc_middles):
				curr_ax = self._axs[i]
				if self.checkbox_coupled.isChecked() == True or self._lastclickedplot == i:
					new_range = [calc_middle-self._range/2, calc_middle+self._range/2]
				else:
					new_range = curr_ax.get_xlim()
				if np.isnan(new_range[0]) or np.isnan(new_range[1]):
					new_range = [0,10000]
				elif new_range[0] == new_range[1]:
					new_range = [new_range[1]-1, new_range[1]+1]
				x_ranges[i] = new_range
				if (self.main_window.automatic_draw_action.isChecked() or dict.get("Force_Draw")) and dict.get("atbu") == None:
					curr_ax.set_xlim(new_range)
			
			timers.append(time.time())
			breakpoint(self, dict, id)

			## y-ranges
			dataframes = []
			files_dicts = []
			with locks["exp_df"]:
				dataframes.append(self._exp_data_df.copy())
				files_dicts.append(self.config["files_experimental"])
			with locks["lin_df"]:
				dataframes.append(self._lin_data_df.copy())
				files_dicts.append(self.config["files_prediction"])
			with locks["ass_df"]:
				dataframes.append(self._catass_df.copy())
				files_dicts.append(self.config["files_assigned"])
				
			for i in axs_to_be_updated:
				breakpoint(self, dict, id)
				datatype = 0
				for dataframe, files in zip(dataframes, files_dicts):
					length = len(dataframe.index)
					if length > 1000:
						x_start = dataframe["x"].searchsorted(x_ranges[i][0], side="left")
						x_stop  = dataframe["x"].searchsorted(x_ranges[i][1], side="right")
						dataframe = dataframe.iloc[x_start:x_stop].copy()

						tmp_len = len(dataframe.index)
						if tmp_len > binning:
							tmp_binning = (x_ranges[i][1]-x_ranges[i][0])/binning
							
							if tmp_binning == 0:
								tmp_binning = tmp_len

							dataframe.loc[:,"binning"] = (dataframe.loc[:,"x"]-x_ranges[i][0])//tmp_binning
							dataframe = dataframe.loc[dataframe.sort_values(["y"]).drop_duplicates("binning", keep="last").sort_values(["x"]).index]
					
					visible_files = {file for file in files.keys() if not files[file].get("hidden", False)}
					# Special Case Hide/Show catalogue files
					if datatype == 2 and self.config["plot_hidecatalogue"] == False:
						visible_files.add("__catalogue__")
					
					if len(visible_files) != len(files) + (datatype == 2):
						dataframe.query("filename in @visible_files", inplace=True)
					xs = dataframe["x"].to_numpy()
					ys = dataframe["y"].to_numpy()
					
					if datatype == 0:
						y_ranges_exp[i] = [dataframe["y"].min(),dataframe["y"].max()]
						segs = (((xs[i],ys[i]),(xs[i+1],ys[1+i])) for i in range(len(xs)-1))
						colors = self.create_colors(files, dataframe)
						exp_coll = matplotlib.collections.LineCollection(segs, colors=colors)
						if "exp" in self._lines[i]:
							self._lines[i]["exp"].remove()
						self._lines[i]["exp"] = self._axs[i].add_collection(exp_coll)
					elif datatype == 1:
						y_ranges_lin[i] = [dataframe["y"].min(),dataframe["y"].max()]
						segs = (((xs[i],0),(xs[i],ys[i])) for i in range(len(xs)))
						colors = self.create_colors(files, dataframe, i)
						lin_coll = matplotlib.collections.LineCollection(segs, colors=colors)
						if "lin" in self._lines[i]:
							self._lines[i]["lin"].remove()
						self._lines[i]["lin"] = self._taxs[i].add_collection(lin_coll)
					elif datatype == 2:
						tuples = list(zip(xs,ys))
						tuples = tuples if len(tuples)!=0 else [[None,None]]
						self._lines[i]["ass"].set_offsets(tuples)
					datatype += 1

			timers.append(time.time())
			breakpoint(self, dict, id)

			if dict.get("atbu") == None:
				if self.config["plot_yrange_exp"] == 0: #highest value of complete data
					min, max = self._yrange
					zoom = self.expScaleZoomInput.value()
					min = min*np.exp(-zoom)
					max = max*np.exp(-zoom)
					tmp_range = (min, max)
					y_ranges_exp = {i:tmp_range for i in axs_to_be_updated}
				elif self.config["plot_yrange_exp"] == 1: #highest value of current data
					y_range = [np.nanmin([y_ranges_exp[x][0] for x in y_ranges_exp]), np.nanmax([y_ranges_exp[x][1] for x in y_ranges_exp])]
					min, max = y_range
					zoom = self.expScaleZoomInput.value()
					min = min*np.exp(-zoom)
					max = max*np.exp(-zoom)
					tmp_range = (min, max)
					y_ranges_exp = {i:tmp_range for i in axs_to_be_updated}
				elif self.config["plot_yrange_exp"] == 2: #make each plot fit the full height
					pass

				if self.config["plot_yrange_lin"] == 0: #highest value of complete data
					min, max = self._tyrange
					zoom = self.linScaleZoomInput.value()
					min = min*np.exp(-zoom)
					max = max*np.exp(-zoom)
					tmp_range = (min, max)
					y_ranges_lin = {i:tmp_range for i in axs_to_be_updated}
				elif self.config["plot_yrange_lin"] == 1: #highest value of current data
					y_range = [np.nanmin([y_ranges_lin[x][0] for x in y_ranges_lin]), np.nanmax([y_ranges_lin[x][1] for x in y_ranges_lin])]
					min, max = y_range
					zoom = self.expScaleZoomInput.value()
					min = min*np.exp(-zoom)
					max = max*np.exp(-zoom)
					tmp_range = (min, max)
					y_ranges_lin = {i:tmp_range for i in axs_to_be_updated}
				elif self.config["plot_yrange_lin"] == 2: #make each plot fit the full height
					pass


				timers.append(time.time())
				breakpoint(self, dict, id)

				## y-ranges
				## set y-ranges
				for i in axs_to_be_updated:
					tmp_range = y_ranges_exp[i]
					if np.isnan(tmp_range[0]) or np.isnan(tmp_range[1]):
						tmp_range = [-1,+1]
					elif tmp_range[0] == tmp_range[1]:
						tmp_range = [tmp_range[0]-1, tmp_range[0]+1]
					else:
						tmp_range = [tmp_range[0]-self.config["plot_ymargin"]*(tmp_range[1]-tmp_range[0]), tmp_range[1]+self.config["plot_ymargin"]*(tmp_range[1]-tmp_range[0])]
						y_ranges_exp[i] = tmp_range
					self._axs[i].set_ylim(tmp_range)

				for i in axs_to_be_updated:
					tmp_range = y_ranges_lin[i]
					if np.isnan(tmp_range[0]) or np.isnan(tmp_range[1]):
						tmp_range = [-1,+1]
					else:
						if y_ranges_exp[i][0] < 0:
							upper = tmp_range[1]+self.config["plot_ymargin"]*(tmp_range[1]-0)
							lower = upper*y_ranges_exp[i][0]/y_ranges_exp[i][1]
						else:
							upper = tmp_range[1]+self.config["plot_ymargin"]*(tmp_range[1]-0)
							lower = 0
						tmp_range = [lower, upper]
					self._taxs[i].set_ylim(tmp_range)

			timers.append(time.time())
			breakpoint(self, dict, id)


			## Tight Layout Experimental
			if self.config["plot_tight_layout"] == True:
				self.fig.tight_layout()
				if self.checkbox_xaxis.isChecked()==False:
					self.fig.subplots_adjust(hspace=0)

			if self.main_window.automatic_draw_action.isChecked() or dict.get("Force_Draw"):
				self.plot_canvas.draw()
			timers.append(time.time())
		
		except CustomError as E:
			self._queue_dict.update({k:dict.get(k, False) for k in ["Force_Draw"]})
		
		finally:
			if self.config["flag_time_plots"] == True:
				deltas = [f"Measured the runtime of 'update plots':"]
				for i in range(len(timers)-1):
					deltas.append(f"Delta {i}: {timers[i+1]-timers[i]:.4f}")
				deltas.append(f"Resulting total runtime: {timers[-1]-timers[0]:.4f}")
				message = '\n    '.join(deltas)
				self.notification(message)


	##
	## Helper for Setting Values
	##
	@synchronized(locks["cat_df"])
	@synchronized(locks["ass_df"])
	def assigned_dataframe(self):
		xs = []
		ys = []
		QNs = []
		tmp_catalogue = self._catalogue.rename(columns={"x": "dump1", "y": "dump2", "freq_predicted": "x", "intens_predicted": "y"})
		tmp_catalogue["filename"] = "__catalogue__"
		tmp_assigned = self._ass_data_df.copy()
		tmp_assigned["y"] = 0
		
		df_join = pd.concat((tmp_catalogue, tmp_assigned), join="inner", ignore_index=True)
		
		qns_labels = ["qnu1", "qnu2", "qnu3", "qnu4", "qnu5", "qnu6", "qnl1", "qnl2", "qnl3", "qnl4", "qnl5", "qnl6"]
		visible = [not x.isHidden() for x in self.QNcbs]
		qns_labels = [x for i,x in enumerate(qns_labels) if visible[i%len(visible)]]
		
		df_join.drop_duplicates(subset=qns_labels, keep="first")
		self._catass_df = df_join.sort_values("x")		
		
	def offsets(self, indices):
		ann_dict = self.config["plot_annotations_dict"]
		len_axs = len(self._axs)
		annotate = self.checkbox_annotate.isChecked()
		series = self._series_method
		
		inds = [len_axs-i-1 for i in indices]
		
		if self._series_method == 0:
			J0 = np.float64(self.input_J0.text())
			try:
				expr = self.inputExpr.text()
				vals = [eval(expr,{"J":ind, "J_0":J0}) for ind in inds]
			except Exception as E:
				vals = [0]*len(indices)
				self.notification("Could not evaluate your expression.")
		
		elif self._series_method == 1:
			vals = []
			i_tuples = self.get_QNs_by_plots(indices)
			files = self.config["files_prediction"]
			visible_files = {file for file in files.keys() if not files[file].get("hidden", False)}
			df = self._lin_data_df.copy()
			if len(visible_files) != len(files):
				df = df.query("filename in @visible_files")
			for tmp_tuple in i_tuples:
				lower, upper, visible, query_string = tmp_tuple
				data = df.query(query_string)["x"].to_numpy()
				val = data[0] if len(data) !=0 else 0
				if self.zeroAssigned.isChecked():
					tmp_val = self._catalogue.query(query_string)["x"].to_numpy()
					if len(tmp_val) != 0:
						val = tmp_val[0]
				vals.append(val)
			
				
		elif self._series_method == 2:
			vals = []
			ind_min = int(self.inputIndexMin.text())
			for ind in inds:
				if ind_min+ind < len(self._frequencies_list):
					val = self._frequencies_list[ind_min+ind]
				else:
					val = 0
				vals.append(val)
		else:
			vals = [0]*len(indices)
			
		for i in indices:
			if i in self._annotations:
				self._annotations[i].set_visible(False)
		if self.checkbox_annotate.isChecked():
			files = self.config["files_assigned"]
			visible_files = {file for file in files.keys() if not files[file].get("hidden", False)}
			if self.config["plot_hidecatalogue"] == False:
				visible_files.add("__catalogue__")
			df = self._catass_df.copy()
			if len(visible_files) != len(files) + (self.config["plot_hidecatalogue"] == False):
				df = df.query("filename in @visible_files")
			if self._series_method == 1:
				for i, tmp_tuple in zip(indices, i_tuples):
					lower, upper, visible, query_string = tmp_tuple
					upper_string = ",".join([str(int(upper[i])) for i in range(len(upper)) if visible[i]])
					lower_string = ",".join([str(int(lower[i])) for i in range(len(lower)) if visible[i]])

					if len(df.query(query_string)) > 0:
						self._annotations[i] = self._axs[i].text(**ann_dict, s=r"${} \leftarrow {}$".format(upper_string, lower_string), transform = self._axs[i].transAxes, color=self.config["color_assigned"])
					else:
						self._annotations[i] = self._axs[i].text(**ann_dict, s=r"${} \leftarrow {}$".format(upper_string, lower_string), transform = self._axs[i].transAxes)
			else:
				for i, val in zip(indices, vals):
					self._annotations[i] = self._axs[i].text(**ann_dict, s=f"{val:.2f}", transform = self._axs[i].transAxes)
		return(np.array(vals))
	
	def offset(self, i):
		return(self.offsets([i])[0])

	def change_scale_mode(self,tornot,i):
		if tornot == 0:
			self.config["plot_yrange_exp"] = i
		elif tornot == 1:
			self.config["plot_yrange_lin"] = i

	def create_colors(self, files={}, dataframe = None, i = None):
		if len(files) == 1 and i == None:
			return(files[list(files.keys())[0]].get("color", "#ffffff"))
		if not isinstance(dataframe, pd.DataFrame):
			return([])

		dataframe.reset_index(drop=True, inplace=True)
		tmp_colors = {file:files[file].get("color", "#ffffff") for file in files.keys()}
		filenames = dataframe["filename"]
		colors = filenames.replace(tmp_colors)

		if i != None and self._series_method == 1:
			lower, upper, visible, query_string = self.get_QNs_by_plot(i)
			index = dataframe.query(query_string).index.to_numpy()
			if len(index) != 0:
				for ind in index:
					if ind in range(len(colors)):
						colors[ind] = self.config["color_current"]

		return(colors)

	def on_hover(self, event):
		x = event.xdata
		y = event.ydata
		if self.checkbox_hover.isChecked() and None not in [x, y, event.inaxes]:
			if event.inaxes in self._taxs.values():
				tax = event.inaxes
				tmp_key = None
				for key, value in self._taxs.items():
					if value == tax:
						tmp_key = key
						break
				ax = self._axs.get(tmp_key)
				if ax != None:
					tmp_coord = tax.transData.transform((x,y))
					inv = ax.transData.inverted()
					coord = inv.transform(tmp_coord)
					x, y = coord

			df = self._lin_data_df
			if len(df) == 0:
				text = ""
			else:
				x_column = df.columns.get_loc("x")
				tmp_can_ind  = df["x"].searchsorted(x, side="left")

				tmp_high_val = df.iloc[tmp_can_ind, x_column] if tmp_can_ind != len(df) else np.inf
				tmp_low_val = df.iloc[tmp_can_ind-1, x_column] if tmp_can_ind != 0 else 0

				if (abs(x-tmp_low_val) > self.config["plot_hover_cutoff"] or tmp_low_val == 0) and abs(x-tmp_high_val) > self.config["plot_hover_cutoff"]:
					text = ""

				else:
					if abs(x-tmp_low_val) <= abs(x-tmp_high_val):
						tmp_low = df["x"].searchsorted(tmp_low_val, side="left")
						tmp_transitions = df.iloc[tmp_low:tmp_can_ind]
					else:
						tmp_high = df["x"].searchsorted(tmp_high_val, side="right")
						if tmp_high > len(df):
							tmp_high = len(df)
						tmp_transitions = df.iloc[tmp_can_ind:tmp_high]
					string = []
					visible = [not(x.isHidden()) for x in self.QNcbs]
					labels_upper = ["qnu1", "qnu2", "qnu3", "qnu4", "qnu5", "qnu6"]
					labels_lower = ["qnl1", "qnl2", "qnl3", "qnl4", "qnl5", "qnl6"]
					for index, row in tmp_transitions.iterrows():
						upper = ",".join([str(int(row[labels_upper[i]])) for i in range(len(labels_upper)) if visible[i]])
						lower = ",".join([str(int(row[labels_lower[i]])) for i in range(len(labels_lower)) if visible[i]])
						string.append(f"{upper} ← {lower}")
					text = " || ".join(string)
			self.topLabel.setText(f"@ ({x:.2f}, {y:.2f}){': ' if text != '' else ''}{text}")
		else:
			self.topLabel.setText("")

	def on_range(self, min, max, i):
		self._lastclickedplot = i
		if max-min == 0:
			return
		modifier = self.modifier_present()
		if (modifier == 0 and self.config["flag_reverse_shift"] == False) or (modifier == -1 and self.config["flag_reverse_shift"] == True) or self.main_window.always_fit_action.isChecked() == True:
			self.assign(i, vmin=min, vmax=max)
		elif (modifier == -1 and self.config["flag_reverse_shift"] == False) or (modifier == 0 and self.config["flag_reverse_shift"] == True): #zoom
			zoom_range = max-min
			if zoom_range == 0:
				return
			middle = (max+min)/2 - self.offset(i)
			self.update_plots({"zoom":{"abs":zoom_range},"step":{"abs":middle}})

	@synchronized(locks["cat_df"])
	@synchronized(locks["ser_df"])
	def assign(self, i, vmin=None, vmax=None, freq=None):
		axs_i = i
		ax = self._axs[i]
		if vmin != None and vmax != None:
			querystring = f"{vmin} <= x <= {vmax}"
			tmp_data = self._exp_data_df.query(querystring)
			tmp_x = tmp_data["x"].to_numpy()
			tmp_y = tmp_data["y"].to_numpy()

			if self._fitline != None:
				self._fitline.remove()
				self._fitline = None
			if self._middlepos != None:
				self._middlepos.remove()
				self._middlepos = None
			# del self._fitline
			# del self._middlepos

			if len(tmp_x)<2:
				self.notification("You did select less than two points of your spectrum, this fit will not work.")
				return()
			try:
				color = self.config["color_fitcolor"]
				if self.config["fit_engine"].startswith("2ndderivative"):
					p0 = [(vmin+vmax)/2, np.amax(tmp_y)-np.amin(tmp_y), (vmax-vmin)/5]
					bounds = [	[vmin, 0, 0],
								[vmax, 3*(np.amax(tmp_y)-np.amin(tmp_y)), 10*(vmax-vmin)]
					]
					
					if self.config["fit_engine"].endswith("Gauss"):
						lineshape_func = lambda *x: lineshape("Gauss", 2, *x)
					elif self.config["fit_engine"].endswith("Lorentz"):
						lineshape_func = lambda *x: lineshape("Lorentz", 2, *x)
					elif self.config["fit_engine"].endswith("Voigt"):
						lineshape_func = lambda *x: lineshape("Voigt", 2, *x)
						p0.append(p0[-1])
						bounds[0].append(bounds[0][-1])
						bounds[1].append(bounds[1][-1])
					
					## Hack as it sometimes works only every second time
					try:
						popt, pciv = optimize.curve_fit(lineshape_func, tmp_x, tmp_y, p0=p0, bounds=bounds)
					except Exception as E:
						popt, pciv = optimize.curve_fit(lineshape_func, tmp_x, tmp_y, p0=p0, bounds=bounds)
					zoom_xs = np.linspace(vmin, vmax, 1000)
					zoom_ys = lineshape_func(zoom_xs, *popt)
					ymax = popt[1]
					xmiddle = popt[0]
					self._fitline, = ax.plot(zoom_xs, zoom_ys, color=color, alpha=0.7, linewidth=1)
				elif self.config["fit_engine"] == "polynom":
					try:
						pol_val = np.polyfit(tmp_x, tmp_y, self.config["fit_rankpolynom"])
					except Exception as E:
						pol_val = np.polyfit(tmp_x, tmp_y, self.config["fit_rankpolynom"])
					pol = np.poly1d(pol_val)
					pol_x = np.linspace(vmin, vmax, 1000)
					pol_y = pol(pol_x)
					ymax = np.max(pol_y)
					xmiddle = pol_x[np.argmax(pol_y)]
					self._fitline, = ax.plot(pol_x, pol_y, color="orange")
				elif self.config["fit_engine"] == "polynom_autorank":
					min_res = np.inf
					opt_rank = 0
					for rank in range(min([20,len(tmp_y)])):
						try:
							try:
								pol_val = np.polyfit(tmp_x, tmp_y, rank)
							except Exception as E:
								pol_val = np.polyfit(tmp_x, tmp_y, rank)
							pol = np.poly1d(pol_val)
							res = np.sum((tmp_y - pol(tmp_x))**2)/len(tmp_y)
							if res < min_res:
								min_res = res
								opt_rank = rank
						except Exception as E:
							continue
					pol_val = np.polyfit(tmp_x, tmp_y, opt_rank)
					pol = np.poly1d(pol_val)
					pol_x = np.linspace(vmin, vmax, 1000)
					pol_y = pol(pol_x)
					ymax = np.max(pol_y)
					xmiddle = pol_x[np.argmax(pol_y)]
					self._fitline, = ax.plot(pol_x, pol_y, color="orange")
				self._middlepos = ax.axvline(x=xmiddle, color=color, ls="--", alpha=1, linewidth=1)
			except Exception as E:
				self._middlepos = None
				self._fitline = None
				self.notification(f"The fitting failed with the following error message : {str(E)}")
				return()
		elif freq!= None:
			xmiddle = freq
			ymax = -1

		tmp_dict={
			"degfreed" :		0,
			"elower" :			0,
			"usd" :				0,
			"tag" :				0,
			"qnfmt":			0,
			"freq_predicted" :	-1,
			"intens_predicted" :-1,
		}
		error = 0
		error_obs_calc = -1
		error_field_value = self.error_field.text()
		try:
			error_field_value = float(error_field_value)*self.config["fit_errorfactor"]
		except ValueError:
			error_field_value = 0
		if self._series_method==1:
			lower, upper, visible, query_string = self.get_QNs_by_plot(axs_i)

			cols = self._lin_data_df.columns
			entry = self._lin_data_df.query(query_string)

			if len(entry) > 0:
				for i in ["degfreed", "elower", "usd", "tag", "qnfmt"]:
					tmp_dict[i] = np.float64(entry[i])
				tmp_dict["tag"] = abs(tmp_dict["tag"])
				if error_field_value == -1:
					error_obs_calc = abs(float(entry["x"])-xmiddle)
				tmp_dict["freq_predicted"] = np.float64(entry["x"])
				tmp_dict["intens_predicted"] = np.float64(entry["y"])

		elif self._series_method == 0:
			lower = [-1]*6
			upper = [-1]*6
			if error_field_value == -1:
				error_obs_calc = abs(self.offset(axs_i)-xmiddle)

		elif self._series_method == 2:
			lower = [-1]*6
			upper = [-1]*6
			if error_field_value == -1:
				error_obs_calc = abs(self.offset(axs_i)-xmiddle)

		if error_field_value == -1 and error_obs_calc != -1:
			error = error_obs_calc
		elif error_field_value == -2:
			resp, rc = QInputDialog.getText(self, 'Set error', 'Error:')
			if rc == True:
				try:
					resp = float(resp)
					error = resp
				except ValueError:
					error = 0
			else:
				error = 0
		else:
			if error_field_value < 0:
				error_field_value = 0
			error = error_field_value

		tmp_dict.update({
			"x" : 		xmiddle,
			"y" :		ymax,
			"error" :		error,
		})
		for i in range(6):
			tmp_dict[f"qnu{i+1}"] = upper[i]
			tmp_dict[f"qnl{i+1}"] = lower[i]
		
		if self.main_window.open_windows.get("series_fit_window") == None:
			cat_df = self._catalogue
			table = self.catalogueTable
			table_model = self.catalogueModel
		else:
			tmp_dict = {k:v for k,v in tmp_dict.items() if k in self._seriesfit_columns}
			cat_df = self._seriesfit
			table = self.main_window.open_windows["series_fit_window"].catalogueTable
			table_model = self.main_window.open_windows["series_fit_window"].catalogueModel
		
		if self.main_window.always_assign_action.isChecked():
			cat_df.reset_index(drop=True, inplace = True)
			new_ind = len(cat_df.index)
			cat_df.loc[new_ind] = tmp_dict
			table.selectRow(new_ind)
			table_model.update()
			table.scrollToBottom()
		else:
			self._lastassign = (axs_i, tmp_dict)
		
		self.update_plots({"atbu":[axs_i]})
	
	@synchronized(locks["cat_df"])
	@synchronized(locks["ser_df"])
	def manual_assign(self):
		if self._lastassign != None:
			if self.main_window.open_windows.get("series_fit_window") == None:
				cat_df = self._catalogue
				table = self.catalogueTable
				table_model = self.catalogueModel
			else:
				tmp_dict = {k:v for k,v in tmp_dict.items() if k in self._seriesfit_columns}
				cat_df = self._seriesfit
				table = self.main_window.open_windows["series_fit_window"].catalogueTable
				table_model = self.main_window.open_windows["series_fit_window"].catalogueModel
			axs_i, tmp_dict = self._lastassign
			cat_df.reset_index(drop=True, inplace = True)
			new_ind = len(cat_df.index)
			cat_df.loc[new_ind] = tmp_dict
			table.selectRow(new_ind)
			table_model.update()
			table.scrollToBottom()
			
			self._lastassign = None
			self.update_plots({"atbu":[axs_i]})
	
	def get_QNs_by_plots(self, indices):
		if self._series_method == 1:
			qnu = [float(x.text()) if x.text()!= "" else 0 for x in self.QNUs]
			qnl = [float(x.text()) if x.text()!= "" else 0 for x in self.QNLs]
			incr = [1 if x.isChecked() else 0 for x in self.QNcbs]
			diff = [1,1,1,1,1,1]
			visible = [not x.isHidden() for x in self.QNcbs]
			len_axs = len(self._axs)
			
			labels_upper = ["qnu1", "qnu2", "qnu3", "qnu4", "qnu5", "qnu6"]
			labels_lower = ["qnl1", "qnl2", "qnl3", "qnl4", "qnl5", "qnl6"]
			
			result = []
			for tmp_ind in indices:
				ind = len_axs-tmp_ind-1
				lower = [qnl[i]+ind*diff[i]*incr[i] if visible[i] else -1 for i in range(len(qnu))]
				upper = [qnu[i]+ind*diff[i]*incr[i] if visible[i] else -1 for i in range(len(qnu))]
				conditions_upper = [f"({labels_upper[i]} == {upper[i]})" for i in range(len(labels_upper)) if visible[i]]
				conditions_lower = [f"({labels_lower[i]} == {lower[i]})" for i in range(len(labels_lower)) if visible[i]]
				query_string = " & ".join(conditions_upper+conditions_lower)
				
				result.append([lower, upper, visible, query_string])
			return(result)
		else:
			return([[([-1]*6, [-1]*6, [False]*6, None)]]*len(indices))
			
	def get_QNs_by_plot(self, i):
		return(self.get_QNs_by_plots([i])[0])

	def alter_QNs(self, index=None, absolute=None):
		if index != None:
			self.config["series_qns"] += index
		elif absolute != None:
			self.config["series_qns"] = absolute
		if self.config["series_qns"] > 6:
			self.config["series_qns"] = 6
		elif self.config["series_qns"] < 1:
			self.config["series_qns"] = 1
		for i in range(self.config["series_qns"]):
			self.QNUs[i].setHidden(False)
			self.QNLs[i].setHidden(False)
			self.QNcbs[i].setHidden(False)
			self.labels[i].setHidden(False)

		for i in range(self.config["series_qns"], len(self.labels)):
			self.QNUs[i].setHidden(True)
			self.QNLs[i].setHidden(True)
			self.QNcbs[i].setHidden(True)
			self.labels[i].setHidden(True)

	@working_d
	@synchronized(locks["axs"])
	def alter_plots(self, index=None, absolute=None, update=True):
		if index == -1 and len(self._axs)==0:
			return
		if index == +1:
			n = len(self._axs)
			if n != 0:
				for i in range(n):
					self._axs[i].change_geometry(n+1, 1, i+1)
					self._taxs[i].change_geometry(n+1, 1, i+1)
			new_ax = self.fig.add_subplot(n+1, 1, n+1)
			# Fill with empty plots, set data in update_plots
			xs = ys = []
			colors = []
			segs = (((xs[i],0),(xs[i],ys[i])) for i in range(len(xs)))

			i = n
			self._axs[i] = new_ax
			self._taxs[i] = new_ax.twinx()
			self._lines[i] = {}
			exp_coll = matplotlib.collections.LineCollection(segs, colors=self.config["color_experimental"])
			self._lines[i]["exp"] = self._axs[i].add_collection(exp_coll)

			lin_coll = matplotlib.collections.LineCollection(segs, colors=colors)
			self._lines[i]["lin"] = self._taxs[i].add_collection(lin_coll)

			self._lines[i]["ass"] = self._taxs[i].scatter(xs, ys, color=self.config["color_assigned"], marker="*")
			rectprops = dict(facecolor='blue', alpha=0.5)
			self._spanSelectors[i] = (matplotlib.widgets.SpanSelector(new_ax, lambda vmax, vmin, index=i:self.on_range(vmax, vmin, index), 'horizontal',rectprops=rectprops, useblit=True))
			self._spanSelectors[i] = (matplotlib.widgets.SpanSelector(self._taxs[i], lambda vmax, vmin, index=i:self.on_range(vmax, vmin, index), 'horizontal',rectprops=rectprops, useblit=True))
		elif index == -1:
			self.fig.delaxes(self._axs[len(self._axs)-1])
			self.fig.delaxes(self._taxs[len(self._taxs)-1])
			del self._axs[len(self._axs)-1]
			del self._taxs[len(self._taxs)-1]
		elif absolute != None:
			for i in range(len(self._axs)):
				self.fig.delaxes(self._axs[i])
				self.fig.delaxes(self._taxs[i])
				del self._axs[i]
				del self._taxs[i]
			self._taxs = {}
			self._annotations = {}
			self._spanSelectors = {}
			self._lines = {}
			absolute = int(absolute)
			if absolute < 1:
				absolute = 1
			self._axs = self.fig.subplots(absolute, gridspec_kw = self.config["plot_matplotlibkwargs"])
			if absolute == 1:
				self._axs = [self._axs]
			self._axs = {i:self._axs[i] for i in range(len(self._axs))}
			xs = ys = []
			segs = (((xs[i],ys[i]),(xs[i+1],ys[1+i])) for i in range(len(xs)-1))
			rectprops = dict(facecolor='blue', alpha=0.5)
			for i in range(absolute):
				curr_ax = self._axs[i]
				self._lines[i] = {}
				self._taxs[i] = curr_ax.twinx()
				exp_coll = matplotlib.collections.LineCollection(segs, colors=self.config["color_experimental"])
				self._lines[i]["exp"] = self._axs[i].add_collection(exp_coll)
				lin_coll = matplotlib.collections.LineCollection(segs)
				self._lines[i]["lin"] = self._taxs[i].add_collection(lin_coll)
				self._lines[i]["ass"] = self._taxs[i].scatter(xs, ys, color=self.config["color_assigned"], marker="*")
				self._spanSelectors[i] = (matplotlib.widgets.SpanSelector(curr_ax, lambda vmax, vmin, index=i:self.on_range(vmax, vmin, index), 'horizontal',rectprops=rectprops, useblit=True))
				self._spanSelectors[i] = (matplotlib.widgets.SpanSelector(self._taxs[i], lambda vmax, vmin, index=i:self.on_range(vmax, vmin, index), 'horizontal',rectprops=rectprops, useblit=True))
		if update == True:
			self.update_plots()

	def QNs_incdec(self, dir):
		upper = (self.QNU1, self.QNU2, self.QNU3)
		lower = (self.QNL1, self.QNL2, self.QNL3)
		checkboxes = (self.QNcb1, self.QNcb2, self.QNcb3)

		for i in range(len(upper)):
			if checkboxes[i].isChecked() == True:
				upper[i].setValue(int(upper[i].text())+dir)
				lower[i].setValue(int(lower[i].text())+dir)
		self.update_plots()

	def series_method(self, x):
		self._series_method = x
		if x == 1 or x == 2: #follow quantum numbers, therfore no offset
			self._x_0 = np.float64(0)

	@synchronized(locks["cat_df"])
	def catalogueTableDelete(self):
		selected = [index.row() for index in self.catalogueTable.selectionModel().selectedRows()]
		for index in sorted(selected, reverse = True):
			self._catalogue.drop(index, inplace =True)
			self._catalogue.reset_index(inplace = True, drop = True)
		self.catalogueModel.update()
		self.update_plots()

	def load_frequencies_list(self, temp = False):
		frequencies = []
		if type(temp) == list:
			frequencies = temp
		elif temp == False:
			fnames = QFileDialog.getOpenFileNames(self, 'Open Frequencies List(s)',)[0]
			if len(fnames)==0:
				return
			for fname in fnames:
				with open(fname, "r") as file:
					for line in file:
						if line.strip() == "" or line.startswith("#"):
							continue
						else:
							tmp_frequencies = re.split('; |, |\s', line)
							for frequency in tmp_frequencies:
								try:
									float_freq = float(frequency)
									frequencies.append(float_freq)
								except ValueError:
									self.notification(f"Could not convert the string '{frequency}' to a numerical value.")
		else:
			line, ok = QInputDialog().getMultiLineText(self, "Specify Custom List","Write list here (delimiters are all whitespace characters, comma and semicolon):")
			if ok and line:
				tmp_frequencies = re.split('; |, |\s', line)
				for frequency in tmp_frequencies:
					try:
						float_freq = float(frequency)
						frequencies.append(float_freq)
					except ValueError:
						self.notification(f"Could not convert the string '{frequency}' to a numerical value.")
		frequencies.sort()
		self._frequencies_list = frequencies

		table = self.frequenciesTable
		table.setRowCount(0)
		i=0
		for frequency in frequencies:
			currRowCount = table.rowCount()
			table.insertRow(currRowCount)

			table.setItem(currRowCount, 0, QTableWidgetItem(f"{i}"))
			table.setItem(currRowCount, 1, QTableWidgetItem(f"{frequency:.4f}"))
			i+=1

	def notification(self, text = None):
		time_str = time.strftime("%H:%M", time.localtime())
		if text == None:
			output = ""
		else:
			output = f"{time_str}: {text}"
		if self.config["flag_debug"] == True:
			print(output)
		self.signalclass.writeLog.emit(output)


	##
	## Creating different Boxes
	##
	def createPlotBox(self):
		## Create Plot Canvas
		self.plotBox = QGroupBox("Plot")
		self.plotBox.setMinimumHeight(200)
		self.plotBox.setMinimumWidth(200)

		self.fig = figure.Figure(figsize=(2.5, 1.5), dpi = self.config["plot_dpi"])
		self.plot_canvas = FigureCanvas(self.fig)
		self.fig.canvas.mpl_connect("motion_notify_event", self.on_hover)

		## Toolbar for debugging
		self.mpl_toolbar = NavigationToolbar2QT(self.plot_canvas, self.plotBox)
		self.mpl_toolbar.setHidden(True)

		## Set Layout
		layout = QVBoxLayout()
		layout.addWidget(self.plot_canvas)
		layout.addWidget(self.mpl_toolbar)
		self.plotBox.setLayout(layout)

		## Initialize different variables
		self._x_0 = np.float64(0)
		self._range = np.float64(10)
		self._axs = self.fig.subplots(0, gridspec_kw = self.config["plot_matplotlibkwargs"])
		self._axs = {i:self._axs[i] for i in range(len(self._axs))}
		self._taxs = {}
		self._annotations = {}
		self._spanSelectors = {}
		self._lines = {}
		self._yrange = [0,1]
		self._tyrange = [0,1]

		## Create and Set up Experimental Dataframe
		##
		self._exp_data_df_columns = ["x", "y", "filename"]
		self._exp_data_df = pd.DataFrame(columns = self._exp_data_df_columns, dtype=np.float64)

		## Create and Set up Prediction Dataframe
		##
		self._lin_data_df_columns = ["x", "error", "y", "degfreed", "elower", "usd", "tag", "qnfmt", 'qnu1', 'qnu2', 'qnu3', 'qnu4', 'qnu5', 'qnu6', 'qnl1', 'qnl2', 'qnl3', 'qnl4', 'qnl5', 'qnl6', "filename"]
		self._lin_data_df_dtypes = [np.float64, np.float64, np.float64, np.int16, np.float64, np.int16, np.int32, np.int16]+[np.int16]*12+[str]
		dtypes_dict = {self._lin_data_df_columns[i]:self._lin_data_df_dtypes[i] for i in range(len(self._lin_data_df_columns))}
		self._lin_data_df = pd.DataFrame(columns = self._lin_data_df_columns)
		self._lin_data_df = self._lin_data_df.astype(dtypes_dict)
	
		## Create and Set up Assigned Dataframe
		##
		self._ass_data_df_columns = ['qnu1', 'qnu2', 'qnu3', 'qnu4', 'qnu5', 'qnu6', 'qnl1', 'qnl2', 'qnl3', 'qnl4', 'qnl5', 'qnl6', "x", "error", "weight", "filename"]
		self._ass_data_df_dtypes = [np.int16]*12 + [np.float64, np.float64, np.float64]+[str]
		dtypes_dict = {self._ass_data_df_columns[i]:self._ass_data_df_dtypes[i] for i in range(len(self._ass_data_df_columns))}
		self._ass_data_df = pd.DataFrame(columns = self._ass_data_df_columns)
		self._ass_data_df = self._ass_data_df.astype(dtypes_dict)

		
		## Create and Set up Catalogue Dataframe
		##
		self._catalogue_columns = ["x", "error", "y", "degfreed", "elower", "usd", "tag", "qnfmt", "qnu1", "qnu2", "qnu3", "qnu4", "qnu5", "qnu6","qnl1", "qnl2", "qnl3", "qnl4", "qnl5", "qnl6", "freq_predicted", "intens_predicted"]
		self._catalogue = pd.DataFrame(columns = self._catalogue_columns)
		self._lastassign = None

		## Create and Set up Cat + Ass Dataframe
		##
		self._catass_columns = list(set(self._catalogue_columns).intersection(set(self._ass_data_df_columns)))+["y", "filename"]
		self._catass_df = pd.DataFrame(columns = self._catass_columns)

		## Create and Set up Series Fit Dataframe
		##
		self._seriesfit_columns = ["x", "qnu1", "qnu2", "qnu3", "qnu4", "qnu5", "qnu6","qnl1", "qnl2", "qnl3", "qnl4", "qnl5", "qnl6"]
		self._seriesfit = pd.DataFrame(columns = self._seriesfit_columns)

		self._fitline = None
		self._middlepos = None
		self._lastclickedplot = 0
		self._assigned_array = []
		self._frequencies_list = []
		self._flag_plots_spaced = True
		self._queue_dict = {}
		self._alter_plots_queue = None

	def createTopLeftGroupBox(self):
		self.topLeftGroupBox = QGroupBox("Configure Plots")
		vboxs = [QVBoxLayout() for i in range(3)]


		vboxs[0].addWidget(QLabel("Exp:"))
		
		self.expScaleZoomInput = QDoubleSpinBox()
		self.expScaleZoomInput.setMaximum(100000)
		self.expScaleZoomInput.setMinimum(0)
		self.expScaleZoomInput.setSingleStep(0.1)
		self.expScaleZoomInput.valueChanged.connect(lambda a: self.update_plots())
		self.expScaleZoomInput.setMaximumWidth(50)
		vboxs[0].addWidget(self.expScaleZoomInput)
		
		self.radio_yexp = QButtonGroup(vboxs[0])
		radioButtons = [QRadioButton(label) for label in ["Full Scale", "Dynamic Scale", "Per Plot Scale"]]
		for i, radioButton in enumerate(radioButtons):
			radioButton.toggled.connect((lambda x: self.change_scale_mode(0, i)))
			self.radio_yexp.addButton(radioButton)
			vboxs[0].addWidget(radioButton)
		
		binninglabel = QLabel("Binning: ")
		self.binning_field = QSpinBox()
		self.binning_field.setMaximumWidth(70)
		self.binning_field.setMaximum(2147483647)
		self.binning_field.setMinimum(1)
		self.binning_field.setSingleStep(100)
		self.binning_field.valueChanged.connect(lambda a: self.update_plots())
		vboxs[0].addWidget(binninglabel)
		vboxs[0].addWidget(self.binning_field)
		
		
		vboxs[1].addWidget(QLabel("Cat:"))
		
		self.linScaleZoomInput = QDoubleSpinBox()
		self.linScaleZoomInput.setMaximum(100000)
		self.linScaleZoomInput.setMinimum(-10000)
		self.linScaleZoomInput.setSingleStep(0.1)
		self.linScaleZoomInput.valueChanged.connect(lambda a: self.update_plots())
		self.linScaleZoomInput.setMaximumWidth(50)
		vboxs[1].addWidget(self.linScaleZoomInput)
		
		self.radio_ylin = QButtonGroup(vboxs[1])
		radioButtons = [QRadioButton(label) for label in ["Full Scale", "Dynamic Scale", "Per Plot Scale"]]
		for i, radioButton in enumerate(radioButtons):
			radioButton.toggled.connect(lambda x: self.change_scale_mode(1, i))
			self.radio_ylin.addButton(radioButton)
			vboxs[1].addWidget(radioButton)
		
		errorlabel = QLabel("Default Error: ")
		self.error_field = QSpinBox()
		self.error_field.setMaximum(10000)
		self.error_field.setMinimum(-2)
		self.error_field.setMaximumWidth(70)
		vboxs[1].addWidget(errorlabel)
		vboxs[1].addWidget(self.error_field)


		self.checkbox_xaxis = QCheckBox("Show x-axis")
		self.checkbox_show_yaxis_exp = QCheckBox("Show y-axis (Exp)")
		self.checkbox_show_yaxis_lin = QCheckBox("Show y-axis (Cat)")
		self.checkbox_coupled = QCheckBox("Plots are coupled")
		self.checkbox_annotate = QCheckBox("Annotate plots")
		self.checkbox_hover = QCheckBox("Show closest transition")


		for widget in [self.checkbox_xaxis, self.checkbox_show_yaxis_exp, self.checkbox_show_yaxis_lin, self.checkbox_coupled, self.checkbox_annotate, self.checkbox_hover]:
			vboxs[2].addWidget(widget)
		
		layout = QHBoxLayout()
		for vbox in vboxs:
			vbox.addStretch(1)
			layout.addLayout(vbox)
		self.topLeftGroupBox.setLayout(layout)

	def createTopRightGroupBox(self):
		self.topRightGroupBox = QGroupBox("Define Series")

		self.tabs = QTabWidget()
		self.tabs.currentChanged.connect(lambda x: self.series_method(x))
		
		tabs = [QWidget() for x in range(3)]
		self.tab2 = tabs[1]
		for i, label in enumerate(("Expression", "Transition", "List")):
			self.tabs.addTab(tabs[i], label)

		# Tab 1 Expression
		tmplayout = QHBoxLayout()
		vboxs = [QVBoxLayout() for i in range(2)]

		labelExpr = QLabel("Expression:")
		labelExpr.setMinimumHeight(10)
		self.inputExpr = QLineEdit()
		self.inputExpr.setMinimumWidth(300)
		self.buttonExpr = QPushButton("Apply")
		self.buttonExpr.setMaximumWidth(150)
		self.buttonExpr.clicked.connect(lambda x: self.update_plots({"apply_series" : True, "step":{"abs":0}}))
		
		labelJ_0 = QLabel("J_0")
		labelJ_0.setMinimumHeight(10)
		self.input_J0 = QSpinBox()
		self.input_J0.setMaximum(1000)
		self.input_J0.valueChanged.connect(lambda x: self.update_plots())
		self.input_J0.setMaximumWidth(50)

		for widget in (labelExpr, self.inputExpr, self.buttonExpr):
			vboxs[0].addWidget(widget)
		for widget in (labelJ_0, self.input_J0, self.buttonExpr):
			vboxs[1].addWidget(widget)
		
		vboxs[1].addStretch(1)

		tmplayout.addLayout(vboxs[0])
		tmplayout.addLayout(vboxs[1])
		tmplayout.addStretch(1)
		tabs[0].setLayout(tmplayout)


		# Tab 2 QNs
		tmptab2 = QHBoxLayout()
		tmplayout = QGridLayout()
		self.QNU1 = QSpinBox()
		self.QNU2 = QSpinBox()
		self.QNU3 = QSpinBox()
		self.QNU4 = QSpinBox()
		self.QNU5 = QSpinBox()
		self.QNU6 = QSpinBox()

		self.QNUs = [self.QNU1, self.QNU2, self.QNU3, self.QNU4, self.QNU5, self.QNU6]

		self.QNL1 = QSpinBox()
		self.QNL2 = QSpinBox()
		self.QNL3 = QSpinBox()
		self.QNL4 = QSpinBox()
		self.QNL5 = QSpinBox()
		self.QNL6 = QSpinBox()

		self.QNLs = [self.QNL1, self.QNL2, self.QNL3, self.QNL4, self.QNL5, self.QNL6]

		self.QNcb1 = QCheckBox("Incr")
		self.QNcb2 = QCheckBox("Incr")
		self.QNcb3 = QCheckBox("Incr")
		self.QNcb4 = QCheckBox("Incr")
		self.QNcb5 = QCheckBox("Incr")
		self.QNcb6 = QCheckBox("Incr")

		self.QNcbs = [self.QNcb1, self.QNcb2, self.QNcb3, self.QNcb4, self.QNcb5, self.QNcb6]

		self.IncQN = QPushButton("Incr")
		self.DecQN = QPushButton("Decr")
		self.IncQN.clicked.connect(lambda x: self.QNs_incdec(+1))
		self.DecQN.clicked.connect(lambda x: self.QNs_incdec(-1))

		self.IncNumberQNs = QPushButton("+")
		self.DecNumberQNs = QPushButton("-")
		self.IncNumberQNs.clicked.connect(lambda x: self.alter_QNs(+1))
		self.DecNumberQNs.clicked.connect(lambda x: self.alter_QNs(-1))

		self.buttonQN = QPushButton("Apply")
		self.buttonQN.clicked.connect(lambda x: self.update_plots({"step":{"abs":0}, "apply_series": True}))

		self.zeroAssigned = QCheckBox("Zero Assigned")

		self.labels = (QLabel("QN 1"), QLabel("QN 2"), QLabel("QN 3"), QLabel("QN 4"), QLabel("QN 5"), QLabel("QN 6"))
		for i in range(len(self.labels)):
			tmplayout.addWidget(self.labels[i], 1, i+1)

		for i in range(len(self.QNUs)):
			tmplayout.addWidget(self.QNUs[i], 2, i+1)
			self.QNUs[i].setMaximumWidth(50)
			self.QNUs[i].setMaximum(10000)

		for i in range(len(self.QNLs)):
			tmplayout.addWidget(self.QNLs[i], 3, i+1)
			self.QNLs[i].setMaximumWidth(50)
			self.QNLs[i].setMaximum(10000)

		for i in range(len(self.QNcbs)):
			tmplayout.addWidget(self.QNcbs[i], 4, i+1)
			self.QNcbs[i].setMaximumWidth(50)

		tmplayout.addWidget(self.IncQN, 2, 7, 1, 2)
		tmplayout.addWidget(self.DecQN, 3, 7, 1, 2)

		tmplayout.addWidget(self.IncNumberQNs, 1, 7)
		tmplayout.addWidget(self.DecNumberQNs, 1, 8)

		tmplayout.addWidget(self.buttonQN, 5, 1, 1, 2)
		tmplayout.addWidget(self.zeroAssigned, 5, 7, 1, 2)

		tmptab2.addLayout(tmplayout)
		tmptab2.addStretch(1)
		tabs[1].setLayout(tmptab2)

		# Tab 3 List
		tmplayout = QGridLayout()

		buttonOpenList = QPushButton("Open List")
		buttonOpenList.clicked.connect(lambda x: self.load_frequencies_list())

		buttonWriteList = QPushButton("Write List")
		buttonWriteList.clicked.connect(lambda x: self.load_frequencies_list(temp=True))

		labelMinInd = QLabel("Min Index: ")

		self.inputIndexMin = QSpinBox()
		self.inputIndexMin.setMaximum(1000000000)
		self.inputIndexMin.setMaximumWidth(50)
		self.inputIndexMin.textChanged.connect(lambda a: self.update_plots())

		listApply = QPushButton("Apply")
		listApply.clicked.connect(lambda x: self.update_plots({"step":{"abs":0}, "apply_series": True}))

		hbox = QHBoxLayout()
		hbox.addWidget(buttonOpenList)
		hbox.addWidget(buttonWriteList)
		hbox.addStretch(1)
		hbox.addWidget(labelMinInd)
		hbox.addWidget(self.inputIndexMin)
		hbox.addWidget(listApply)

		tmplayout.addLayout(hbox, 0, 0, 1, 3)

		self.frequenciesTable = QTableWidget()
		self.frequenciesTable.setMaximumHeight(80)
		self.frequenciesTable.setRowCount(0)
		self.frequenciesTable.setColumnCount(2)
		self.frequenciesTable.move(0,0)
		self.frequenciesTable.setHorizontalHeaderLabels(["#", "Frequency"])
		tmplayout.addWidget(self.frequenciesTable, 1, 0, 3, 3)

		tabs[2].setLayout(tmplayout)

		toplayout = QHBoxLayout()

		noplabel = QLabel("Plots")
		self.series_nop_input = QSpinBox()
		self.series_nop_input.setMaximum(500)
		self.series_nop_input.setMinimum(1)
		self.series_nop_input.setMaximumWidth(100)

		widthlabel = QLabel("Width")
		self.series_width_input = QSpinBox()
		self.series_width_input.setMaximum(100000000)
		self.series_width_input.setSingleStep(20)
		self.series_width_input.setMaximumWidth(50)

		self.series_ignore_checkbox = QCheckBox("Ignore on apply")

		toplayout.addWidget(noplabel)
		toplayout.addWidget(self.series_nop_input)
		toplayout.addWidget(widthlabel)
		toplayout.addWidget(self.series_width_input)
		toplayout.addStretch(1)
		toplayout.addWidget(self.series_ignore_checkbox)

		layout = QVBoxLayout()
		layout.addLayout(toplayout)
		layout.addWidget(self.tabs)
		layout.addStretch(1)
		self.topRightGroupBox.setLayout(layout)

	def createCatalogueBox(self):
		self.catalogueBox = QGroupBox("Assigned Lines")

		self.catalogueTableDeleteButton = QPushButton("X")
		self.catalogueTableDeleteButton.setDefault(True)
		self.catalogueTableDeleteButton.clicked.connect(lambda x: self.catalogueTableDelete())
		self.catalogueTableDeleteButton.setFixedWidth(20)

		self.catalogueTableResizeButton = QPushButton("Tight")
		self.catalogueTableResizeButton.setDefault(True)
		self.catalogueTableResizeButton.clicked.connect(lambda x: self.catalogueTable.resizeColumnsToContents())
		self.catalogueTableResizeButton.setFixedWidth(50)

		headers = ["Freq", "Err.", "Intens.", "DF", "Ener.", "USD", "Tag", "Q", "U1", "U2", "U3", "U4", "U5", "U6", "L1", "L2", "L3", "L4", "L5", "L6"]

		self.catalogueTable = QTableView()
		self.catalogueModel = CatalogueTableModel(self._catalogue, headers, 2)
		self.catalogueTable.setModel(self.catalogueModel)
	
		layout = QVBoxLayout()
		buttonsBox = QHBoxLayout()
		buttonsBox.addWidget(self.catalogueTableDeleteButton)
		buttonsBox.addStretch(2)
		buttonsBox.addWidget(self.catalogueTableResizeButton)
		layout.addLayout(buttonsBox)
		layout.addWidget(self.catalogueTable)

		self.catalogueBox.setLayout(layout)

	def createLogBox(self):
		self.logBox = QGroupBox("Log")
		layout = QVBoxLayout()

		self.signalclass.writeLog.connect(lambda text: self.writeLog(text))

		self.log_area = QTextEdit()
		self.log_area.setReadOnly(True)
		self.log_area.setMinimumHeight(50)

		layout.addWidget(self.log_area)
		self.logBox.setLayout(layout)

	def writeLog(self, text):
		self.log_area.append(text)
		sb = self.log_area.verticalScrollBar()
		sb.setValue(sb.maximum())
		if self.logBox.isHidden() and self.main_window.always_show_log_action.isChecked():
			self.logBox.setHidden(False)

	def createQuotesBox(self):
		self.quotesBox = QGroupBox("Quote of the day")
		self.quotesBox.setHidden(True)
		layout = QHBoxLayout()

		self.quote = QLabel()
		self.quote.setWordWrap(True)
		self.quote.setAlignment(Qt.AlignCenter)

		layout.addWidget(self.quote)
		self.quotesBox.setLayout(layout)

	def command_line(self):
		text, ok = QInputDialog().getMultiLineText(self, "Console Input","Enter console command:")
		if ok and text:
			message = []
			old_stdout = sys.stdout
			red_output = sys.stdout = StringIO()
			try:
				exec(text)
			except Exception as E:
				message.append(f"Executing the code raised an error: {str(E)}")
			finally:
				sys.stdout = old_stdout

			message.append("\n".join([f">>> {line}" for line in text.split("\n")]))
			message.append(red_output.getvalue())
			self.notification("\n".join(message))
	
	def new_quote(self):
		quotes = json.loads(quotes_str)
		quote = quotes[random.randint(0,len(quotes)-1)]
		self.quote.setText(quote)
		self.quote.setSizePolicy(QSizePolicy.Minimum, QSizePolicy.Minimum)

	##
	## Event Handler
	##
	def modifier_present(self):
		modifiers = QApplication.keyboardModifiers()
		if modifiers == Qt.ShiftModifier:
			return(0) #shift
		elif modifiers == Qt.ControlModifier:
			return(1) #ctrl
		elif modifiers == (Qt.ControlModifier | Qt.ShiftModifier):
			return(2) #ctrl and shift
		else:
			return(-1)

	def wheelEvent(self,event):
		modifier = self.modifier_present()
		steps = event.angleDelta().y() // 120
		if self.plotBox.underMouse():
			steps = 2**(-steps)
			self.update_plots({"zoom":{"factor":steps}})
		elif self.tab2.underMouse():
			self.QNs_incdec(steps)

	def keyPressEvent(self, event):
		modifier = self.modifier_present()
		# no modifier pressed
		if modifier == -1:
			if event.key() in [65]:
				self.update_plots({"step":{"dir":-1}})
			elif event.key() in [68]:
				self.update_plots({"step":{"dir":+1}})
			elif event.key() in [87]:
				self.update_plots({"zoom":{"factor":1/2}})
			elif event.key() in [83]:
				self.update_plots({"zoom":{"factor":2}})
			elif event.key() in [43]:
				self.alter_plots(+1)
			elif event.key() in [45]:
				self.alter_plots(-1)
		# shift
		elif modifier == 0:
			if event.key() in [65]:
				self.update_plots({"step":{"percent":-20}})
			elif event.key() in [68]:
				self.update_plots({"step":{"percent":+20}})
			elif event.key() in [87]:
				self.update_plots({"zoom":{"factor":3/4}})
			elif event.key() in [83]:
				self.update_plots({"zoom":{"factor":4/3}})
			elif event.key() in [69]:
				self.QNs_incdec(+1)
			elif event.key() in [81]:
				self.QNs_incdec(-1)
			elif event.key() in [32]:
				self.update_plots({"step":{"abs":0}, "apply_series": True})
			elif event.key() in [75]:
				self.command_line()
			elif event.key() in [61]:
				self.new_quote()
				self.quotesBox.setHidden(not self.quotesBox.isHidden())
		elif modifier == 1:
			pass
		elif modifier == 2:
			pass


##
## Enhanced Window Class
##
class EQWidget(QWidget):
	def __init__(self, id, parent=None):
		self.main_window = mw
		self.main_widget = mw.main_widget
		self.config = self.main_widget.config
		super().__init__(parent)
	
	@synchronized(locks["windows"])
	def closeEvent(self, event):
		self.main_window.open_windows[id] = None

	def keyPressEvent(self, event):
		if event.key() == Qt.Key_Escape:
			self.close()
		self.main_widget.keyPressEvent(event)

##
## Files Window Class
##
class FilesWindow(EQWidget):
	def __init__(self, files, actions, title, df, lock, parent=None):
		super().__init__(parent)
		self.setWindowTitle(title)
		self.files = files
		self.actions = actions
		self.df = df
		self.lock = lock
		self.files_dict = {}
		self.main_widget.signalclass.updateWindows.connect(self.update)
		self.layout = None
		self.update()
		
	def GUI_top(self):
		pass
		
	def GUI_bottom(self):
		pass
	
	def GUI_init(self):
		if self.layout == None:
			self.layout = QVBoxLayout()
			self.layout_actions = QHBoxLayout()
			self.layout_files = QGridLayout()
			self.layout.addLayout(self.layout_actions)
			self.layout.addLayout(self.layout_files)
			self.layout.addStretch(1)
			tmp_button = QPushButton("Load")
			tmp_button.clicked.connect(lambda a: self.FUNC_add())
			self.layout_actions.addWidget(tmp_button)
			tmp_button = QPushButton("Add")
			tmp_button.clicked.connect(lambda a: self.FUNC_add(True))
			self.layout_actions.addWidget(tmp_button)
			self.layout_actions.addStretch(1)
			self.setLayout(self.layout)
			self.setMinimumHeight(200)
			self.setMinimumWidth(600)
		layout = self.layout_files

		self.row_id = 0
		
		self.GUI_top()

		for file in self.files:
			tmp_dict = {	"label" :			QLabel(file),
							"input" :			QLineEdit(),
							"button" :			QPushButton("CP"),
							"button_invert" :	QPushButton("Invert"),
							"button_hide" :		QPushButton("Show" if self.files[file].get("hidden")==True else "Hide"),
							"button_delete" :	QPushButton("X"),
			}
			
			tmp_color = self.files[file].get("color", None)
			
			tmp_dict["input"].setText(tmp_color)
			tmp_dict["input"].textChanged.connect(lambda a,x=file: self.change_color(x, inp = True))
			
			tmp_dict["button"].clicked.connect(lambda a,x=file: self.change_color(x))
			tmp_dict["button"].setStyleSheet(f"background-color: {tmp_color}")
			
			tmp_dict["button_delete"].clicked.connect(lambda a,x=file: self.delete_file(x))
			tmp_dict["button_hide"].clicked.connect(lambda a, x=file: self.hide_file(x))
			tmp_dict["button_invert"].clicked.connect(lambda a, x=file: self.invert_file(x))
			
			col_id = 0
			layout.addWidget(tmp_dict["label"], self.row_id, col_id)
			col_id += 1
			
			for action in self.actions:
				if action == "color":
					layout.addWidget(tmp_dict["input"], self.row_id, col_id)
					col_id += 1
					layout.addWidget(tmp_dict["button"], self.row_id, col_id)
				elif action == "hide":
					layout.addWidget(tmp_dict["button_hide"], self.row_id, col_id)
				elif action == "invert":
					layout.addWidget(tmp_dict["button_invert"], self.row_id, col_id)
				elif action == "delete":
					layout.addWidget(tmp_dict["button_delete"], self.row_id, col_id)
				col_id += 1

			self.files_dict[file] = tmp_dict
			self.row_id += 1
			
		self.GUI_bottom()
		
	def update(self):
		for key, value in self.files_dict.items():
			for widget in value.values():
				widget.setParent(None)
		self.files_dict = {}
		self.GUI_init()
		
	def check_color(self, file, inp = False):
		if inp == False:
			color = QColorDialog.getColor()
			color = color.name()
		else:
			color = self.files_dict[file]["input"].text()
		return(match_color(color))
	
	@synchronized(locks["axs"])
	def hide_file(self, file, inp = False):
		hidden = self.files[file].get("hidden", False)
		hidden = not hidden
		self.files[file]["hidden"] = hidden
		
		if hidden:
			self.files_dict[file]["button_hide"].setText("Show")
		else:
			self.files_dict[file]["button_hide"].setText("Hide")
		
		self.main_widget.update_plots()
	
	@threading_d
	def delete_file(self, file):
		del self.files[file]
		args = {}
		with self.lock:
			if self.df == "lin":
				tmp_df = self.main_widget._lin_data_df
				args["do_QNs"] = False
			elif self.df == "exp":
				tmp_df = self.main_widget._exp_data_df
			else:
				tmp_df = self.main_widget._ass_data_df
			tmp_df.drop(tmp_df[tmp_df["filename"]==file].index, inplace=True)
		self.main_window.load_cat(keep_old = True, add_files = False, **args)

##
## Cat Files Window
##
class LinesWindow(FilesWindow):
	def __init__(self, parent=None):
		title = "Cat Files"
		actions = ["color", "hide", "delete"]
		files = mw.config["files_prediction"]
		df = "lin"
		lock = locks["lin_df"]
		super().__init__(files, actions, title, df, lock, parent)
	
	def FUNC_add(self, keep_old = False):
		self.main_window.load_cat(keep_old = keep_old)
	
	def GUI_top(self):
		for file, title, color in zip(("__starassigned__", "__current__"), ("Assigned Markers", "Current Transition"), (self.config["color_assigned"], self.config["color_current"])):
			tmp_dict = {	"label" :			QLabel(title),
							"input" :			QLineEdit(),
							"button" :			QPushButton("CP"),
			}
			tmp_color = color
			tmp_dict["input"].setText(tmp_color)
			tmp_dict["input"].textChanged.connect(lambda a,x=file: self.change_color(x, inp = True))
			
			tmp_dict["button"].clicked.connect(lambda a,x=file: self.change_color(x))
			tmp_dict["button"].setStyleSheet(f"background-color: {tmp_color}")
			
			col_id = 0
			self.layout_files.addWidget(tmp_dict["label"], self.row_id, col_id)
			col_id += 1
			self.layout_files.addWidget(tmp_dict["input"], self.row_id, col_id)
			col_id += 1
			self.layout_files.addWidget(tmp_dict["button"], self.row_id, col_id)
			col_id += 1
			
			self.files_dict[file] = tmp_dict
			self.row_id += 1
	
	@synchronized(locks["axs"])
	def change_color(self, file, inp = False):
		color = self.check_color(file, inp)
		if color:
			self.files_dict[file]["button"].setStyleSheet(f"background-color: {color}")
			self.files_dict[file]["input"].setText(color)
			if file == "__starassigned__":
				self.config["color_assigned"] = color
				for i in range(len(self.main_widget._lines)):
					self.main_widget._lines[i]["ass"].set_color(color)
			elif file == "__current__":
				self.config["color_current"] = color
			else:
				self.files[file]["color"] = color
			self.main_widget.update_plots()

##
## Spectra Files Window
##
class SpectraWindow(FilesWindow):
	def __init__(self, parent=None):
		title = "Spectrum Files"
		actions = ["color", "hide", "invert", "delete"]
		files = mw.config["files_experimental"]
		df = "exp"
		lock = locks["exp_df"]
		super().__init__(files, actions, title, df, lock, parent)

	def FUNC_add(self, keep_old = False):
		self.main_window.load_exp(keep_old = keep_old)

	def GUI_top(self):
		file = "__defaultexpcolor__"
		tmp_dict = {	"label" :			QLabel("Initial Color"),
						"input" :			QLineEdit(),
						"button" :			QPushButton("CP"),
		}
		tmp_color = self.config.get('color_experimental')
		tmp_dict["input"].setText(tmp_color)
		tmp_dict["input"].textChanged.connect(lambda a,x=file: self.change_color(x, inp = True))
		
		tmp_dict["button"].clicked.connect(lambda a,x=file: self.change_color(x))
		tmp_dict["button"].setStyleSheet(f"background-color: {tmp_color}")
		
		col_id = 0
		self.layout_files.addWidget(tmp_dict["label"], self.row_id, col_id)
		col_id += 1
		self.layout_files.addWidget(tmp_dict["input"], self.row_id, col_id)
		col_id += 1
		self.layout_files.addWidget(tmp_dict["button"], self.row_id, col_id)
		col_id += 1
		
		self.files_dict[file] = tmp_dict
		self.row_id += 1

	def change_color(self, file, inp=False):
		color = self.check_color(file, inp)
		if color:
			self.files_dict[file]["button"].setStyleSheet(f"background-color: {color}")
			self.files_dict[file]["input"].setText(color)
			if file == "__defaultexpcolor__":
				self.config["color_experimental"] = color
			else:
				self.files[file]["color"] = color
			self.main_widget.update_plots()

	@synchronized(locks["exp_df"])
	def invert_file(self, file):
		self.files[file]["invert"] = not self.files[file].get("invert", False)
		tmp_df = self.main_widget._exp_data_df
		tmp_df.loc[tmp_df.filename == file, "y"] = -tmp_df.loc[tmp_df.filename == file, "y"]
		self.main_widget.update_plots()


##
## Assigned Files Window
##
class AssignedWindow(FilesWindow):
	def __init__(self, parent=None):
		title = "Assigned Files"
		actions = ["hide", "delete"]
		files = mw.config["files_assigned"]
		df = "ass"
		lock = locks["ass_df"]
		super().__init__(files, actions, title, df, lock, parent)

	def FUNC_add(self, keep_old = False):
		self.main_window.load_lin(keep_old = keep_old)
	
	def GUI_top(self):
		file = "__catalogue__"
		tmp_dict = {	"label" :			QLabel("Here Assigned Lines"),
						"button_hide" :		QPushButton("Show" if self.config["plot_hidecatalogue"]==True else "Hide"),
		}
		
		tmp_dict["button_hide"].clicked.connect(lambda a,x=file: self.hide_catalogue(x))
		
		col_id = 0
		self.layout_files.addWidget(tmp_dict["label"], self.row_id, col_id)
		col_id += 1
		self.layout_files.addWidget(tmp_dict["button_hide"], self.row_id, col_id)
		self.files_dict[file] = tmp_dict
		self.row_id += 1
	
	@synchronized(locks["axs"])
	def hide_catalogue(self, file, inp = False):
		hidden = not self.config.get("plot_hidecatalogue", False)
		self.config["plot_hidecatalogue"] = hidden
		if hidden:
			self.files_dict[file]["button_hide"].setText("Show")
		else:
			self.files_dict[file]["button_hide"].setText("Hide")
		self.main_widget.update_plots()


##
## Lineshape Window
##
class LineshapeWindow(EQWidget):
	def __init__(self, parent=None):
		super().__init__(parent)
		self.setWindowTitle("Lineshape Window")


		self.main_window = mw
		self.main_widget = mw.main_widget
		self.config = self.main_widget.config
		if self.main_widget._lastclickedplot in self.main_widget._axs.keys():
			self.i = self.main_widget._lastclickedplot
		else:
			self.i = 0
		layout = QGridLayout()

		labelFunction = QLabel("Function: ")
		self.input_func = QComboBox()
		self.input_func.addItem("Gauss")
		self.input_func.addItem("Lorentz")
		self.input_func.addItem("Voigt")

		layout.addWidget(labelFunction, 0, 0)
		layout.addWidget(self.input_func, 0, 1)

		labelWidthGauss = QLabel("Gauss: ")
		self.input_width_gauss = QLineEdit()
		self.input_width_gauss.setMaximumWidth(200)

		layout.addWidget(labelWidthGauss, 0, 2)
		layout.addWidget(self.input_width_gauss, 0, 3)

		labelWidthLorentz = QLabel("Lorentz: ")
		self.input_width_lorentz = QLineEdit()
		self.input_width_lorentz.setMaximumWidth(200)

		layout.addWidget(labelWidthLorentz, 0, 4)
		layout.addWidget(self.input_width_lorentz, 0, 5)

		labelDer = QLabel("Derivative: ")
		self.input_derivative = QSpinBox()
		self.input_derivative.setMaximum(2)
		self.input_derivative.setMaximumWidth(200)

		layout.addWidget(labelDer, 0, 6)
		layout.addWidget(self.input_derivative, 0, 7)

		layout.setColumnStretch(8, 1)

		self.button_apply = QPushButton("Apply")
		layout.addWidget(self.button_apply, 0, 9)

		self.dpi = self.config["plot_dpi"]
		self.fig = figure.Figure(figsize=(2.5, 1.5), dpi = self.dpi)
		self.plot_canvas = FigureCanvas(self.fig)
		layout.addWidget(self.plot_canvas, 1, 0, 4, 12)
		self.setLayout(layout)

		self._ax = self.fig.subplots()
		self._tax = self._ax.twinx()

		self._plot_exp = self._ax.plot([],[], color = self.config["color_experimental"])
		self._plot_lin = self._tax.plot([],[], color = "blue")

		rectprops = dict(facecolor='blue', alpha=0.5)
		self.span_selector1 = matplotlib.widgets.SpanSelector(self._ax, lambda vmax, vmin: self.main_widget.on_range(vmax, vmin, self.i), 'horizontal',rectprops=rectprops)
		self.span_selector2 = matplotlib.widgets.SpanSelector(self._tax, lambda vmax, vmin: self.main_widget.on_range(vmax, vmin, self.i), 'horizontal',rectprops=rectprops)
		
		self.input_func.setCurrentText(self.config["lineshapewindow_lineshape"])
		self.input_derivative.setValue(self.config["lineshapewindow_derivative"])
		self.input_width_gauss.setText(self.config["lineshapewindow_gauss"])
		self.input_width_lorentz.setText(self.config["lineshapewindow_lorentz"])
		
		
		self.plot()
		self.button_apply.clicked.connect(self.plot)

		self.cid = self.main_widget.fig.canvas.mpl_connect("draw_event", lambda a:self.plot())

	def plot(self):
		if self.main_widget._lastclickedplot in self.main_widget._axs.keys():
			self.i = self.main_widget._lastclickedplot
		else:
			self.i = 0
		x_range = self.main_widget._axs[self.i].get_xlim()
		dataframe_exp = self.main_widget._exp_data_df.query(f"{x_range[0]} <= x <= {x_range[1]}")
		exp_range = [dataframe_exp["y"].min(), dataframe_exp["y"].max()]
		dataframe_lin = self.main_widget._lin_data_df.query(f"{x_range[0]} <= x <= {x_range[1]}")
		lin_xs = dataframe_lin["x"].to_numpy()
		lin_ys = dataframe_lin["y"].to_numpy()

		try:
			lorentz = float(self.input_width_lorentz.text())
		except ValueError:
			lorentz = 1
		try:
			gauss = float(self.input_width_gauss.text())
		except ValueError:
			gauss = 1
		try:
			derivative = float(self.input_derivative.text())
		except ValueError:
			derivative = 0
			
		function = self.input_func.currentText()
		xs = np.linspace(x_range[0], x_range[1], 1000)
		args = [function, derivative, xs, 0, 0]
		if function == "Gauss" or function == "Voigt":
			args.append(gauss)
		if function == "Lorentz" or function == "Voigt":
			args.append(lorentz)
		
		yss = []
		for i in range(len(lin_xs)):
			args[3] = lin_xs[i]
			args[4] = lin_ys[i]
			ys = lineshape(*args)
			yss.append(ys)
		if len(yss) > 0:
			ys = np.mean(np.array(yss), axis=0)
		else:
			ys = np.zeros(len(xs))

		lin_range = [ys.min(), ys.max()]

		self._plot_exp[0].set_data(dataframe_exp["x"].to_numpy(), dataframe_exp["y"].to_numpy())
		self._plot_lin[0].set_data(xs, ys)

		margin = self.config["plot_ymargin"]
		if np.isnan(exp_range[0]) or np.isnan(exp_range[1]) or exp_range[0] == exp_range[1]:
			exp_range = [-1,+1]

		if np.isnan(lin_range[0]) or np.isnan(lin_range[1]) or lin_range[0] == lin_range[1]:
			lin_range = [-1,+1]
		else:
			lin_range = [exp_range[0]*lin_range[1]/exp_range[1], lin_range[1]]

		exp_range = [exp_range[0]-margin*(exp_range[1]-exp_range[0]), exp_range[1]+margin*(exp_range[1]-exp_range[0])]
		lin_range = [lin_range[0]-margin*(lin_range[1]-lin_range[0]), lin_range[1]+margin*(lin_range[1]-lin_range[0])]

		self._ax.set_xlim(x_range)
		self._ax.set_ylim(exp_range)
		self._tax.set_ylim(lin_range)
		self.plot_canvas.draw()

	@synchronized(locks["windows"])
	def closeEvent(self, event):
		self.config["lineshapewindow_lineshape"] = self.input_func.currentText()
		self.config["lineshapewindow_derivative"] = self.input_derivative.value()
		self.config["lineshapewindow_gauss"] = self.input_width_gauss.text()
		self.config["lineshapewindow_lorentz"] = self.input_width_lorentz.text()
		self.main_widget.fig.canvas.mpl_disconnect(self.cid)
		super().closeEvent(event)

##
## Blended Peaks Window
##
class BlendedlinesWindow(EQWidget):
	def __init__(self, parent=None):
		super(BlendedlinesWindow, self).__init__(parent)
		self.setWindowTitle("Blended Lines Window")


		self.main_window = mw
		self.main_widget = mw.main_widget
		self.config = self.main_widget.config
		if self.main_widget._lastclickedplot in self.main_widget._axs.keys():
			self.i = self.main_widget._lastclickedplot
		else:
			self.i = 0
		layout = QGridLayout()

		self.peaks = []
		self.x_range = (0,0)
		self.exp_range = (0,0)
		self.polynom = 0
		
		labelFunction = QLabel("Function: ")
		self.input_func = QComboBox()
		self.input_func.addItem("Gauss")
		self.input_func.addItem("Lorentz")
		self.input_func.addItem("Voigt")

		layout.addWidget(labelFunction, 0, 0)
		layout.addWidget(self.input_func, 0, 1)

		labelDer = QLabel("Derivative: ")
		self.input_derivative = QSpinBox()
		self.input_derivative.setMaximum(2)
		self.input_derivative.setMaximumWidth(200)

		layout.addWidget(labelDer, 0, 2)
		layout.addWidget(self.input_derivative, 0, 3)

		self.button_add = QPushButton("Add Peak")
		layout.addWidget(self.button_add, 0, 4)
		self.button_add.clicked.connect(lambda a:self.add_peak())

		self.button_del = QPushButton("Del Peak")
		layout.addWidget(self.button_del, 0, 5)
		self.button_del.clicked.connect(self.del_peak)

		self.button_delall = QPushButton("Del All")
		layout.addWidget(self.button_delall, 0, 6)
		self.button_delall.clicked.connect(lambda a:self.del_peak(-1))

		labelWidth = QLabel("Max FWHM: ")
		self.input_width = QLineEdit()
		self.input_width.setMaximumWidth(50)
		layout.addWidget(labelWidth, 0, 7)
		layout.addWidget(self.input_width, 0, 8)

		labelTrans = QLabel("Transpar.: ")
		self.input_transp = QDoubleSpinBox()
		self.input_transp.setValue(0.10)
		self.input_transp.setMaximumWidth(200)
		self.input_transp.setRange(0, 1)
		self.input_transp.setSingleStep(0.1)
		layout.addWidget(labelTrans, 0, 9)
		layout.addWidget(self.input_transp, 0, 10)

		labelPol = QLabel("Polynom: ")
		self.checkbox_polynom = QCheckBox("Show")
		self.input_polynom = QSpinBox()
		layout.addWidget(labelPol, 0, 11)
		layout.addWidget(self.input_polynom, 0, 12)
		layout.addWidget(self.checkbox_polynom, 0, 13)

		layout.setColumnStretch(14, 1)

		self.button_apply = QPushButton("Update")
		layout.addWidget(self.button_apply, 0, 15)
		self.button_apply.clicked.connect(self.plot)

		self.dpi = self.config["plot_dpi"]
		self.fig = figure.Figure(figsize=(2.5, 1.5), dpi = self.dpi)
		self.plot_canvas = FigureCanvas(self.fig)
		layout.addWidget(self.plot_canvas, 1, 0, 4, 16)
		layout.setRowStretch(4, 1)

		self.table = QTableWidget()
		self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
		self.table.setMinimumHeight(50)
		layout.addWidget(self.table, 5, 0, 1, 16)

		self.setLayout(layout)

		self._ax = self.fig.subplots()

		self._plot_exp = self._ax.plot([],[], color = self.config["color_experimental"])
		self._plot_fit = self._ax.plot([],[], color = "blue")
		self._plots_parts = []

		self.input_func.setCurrentText(self.config["blendedlineswindow_lineshape"])
		self.input_derivative.setValue(int(self.config["blendedlineswindow_derivative"]))
		self.input_width.setText(self.config["blendedlineswindow_maxfwhm"])
		self.input_transp.setValue(self.config["blendedlineswindow_transparency"])
		self.input_polynom.setValue(self.config["blendedlineswindow_polynom"])
		self.checkbox_polynom.setChecked(self.config["blendedlineswindow_showbaseline"])
		
		self.plot()
		self.update()

		self.cid = self.fig.canvas.mpl_connect("button_press_event", lambda event: self.add_peak(event.xdata, event.ydata))

	def add_peak(self, x=None, y=None):
		x_range = self.x_range
		y_range = self.exp_range
		if x == None:
			x = (x_range[0]+x_range[1])/2
		if y == None:
			y = y_range[1]
		try:
			x = np.float64(x)
			y = np.float64(y)
		except ValueError:
			x = (x_range[0]+x_range[1])/2
			y = y_range[1]
		width = (x_range[1]-x_range[0])/10
		self.peaks.append((x,y,width,width))
		self.update()

	def del_peak(self, i=None):
		if i == -1:
			self.peaks = []
		elif i!=None and np.issubdtype(type(i), np.number):
			if i in range(len(self.peaks)):
				del self.peaks[int(i)]
		else:
			if len(self.peaks) != 0:
				self.peaks.pop()
		self.update()

	def f(self, x, *ps, pol=True):
		try:
			derivative = float(self.input_derivative.text())
		except ValueError:
			derivative = 0
		function = self.input_func.currentText()
		signal = np.zeros(len(x))
		if pol != False and self.polynom > 0:
			pol = self.polynom
			pspol = ps[-pol:]
			ps = ps[:-pol]
			signal += np.polyval(pspol, x-sum(self.x_range)/2)
		for i in range(int(len(ps)/self.arg_nr)):
			signal += lineshape(function, derivative, x, *ps[i*self.arg_nr: (i+1)*self.arg_nr])
		return(signal)

	def plot(self):
		if self.main_widget._lastclickedplot in self.main_widget._axs.keys():
			self.i = self.main_widget._lastclickedplot
		else:
			self.i = 0
		self.x_range = x_range = self.main_widget._axs[self.i].get_xlim()
		dataframe_exp = self.main_widget._exp_data_df.query(f"{x_range[0]} <= x <= {x_range[1]}")
		self.exp_range = exp_range = [dataframe_exp["y"].min(), dataframe_exp["y"].max()]
		self.exp_xs = dataframe_exp["x"].to_numpy()
		self.exp_ys = dataframe_exp["y"].to_numpy()

		self._plot_exp[0].set_data(self.exp_xs, self.exp_ys)

		margin = self.config["plot_ymargin"]
		if np.isnan(exp_range[0]) or np.isnan(exp_range[1]) or exp_range[0] == exp_range[1]:
			exp_range = [-1,+1]

		exp_range = [exp_range[0]-margin*(exp_range[1]-exp_range[0]), exp_range[1]+margin*(exp_range[1]-exp_range[0])]

		self._ax.set_xlim(x_range)
		self._ax.set_ylim(exp_range)
		self.update()

	def update(self):
		trans = self.input_transp.value()
		function = self.input_func.currentText()
		self.arg_nr = 4 if function == "Voigt" else 3
		for pl in self._plots_parts:
			pl.remove()
		self._plots_parts = []
		x_range = self.x_range
		if len(self.exp_ys) != 0:
			y_range = [np.amin(self.exp_ys), np.amax(self.exp_ys)]
		else:
			y_range = [0,0]
		xs = np.linspace(x_range[0], x_range[1], 1000)
		
		ys_fit = self.exp_ys
		self.polynom = self.input_polynom.value()
		if len(ys_fit) != 0 and (self.polynom + len(self.peaks)) != 0:
			p0 = []
			bounds = [[],[]]
			try:
				width = np.float64(self.input_width.text())
			except ValueError:
				width = (x_range[1]-x_range[0])/2
			amp = 4*(np.amax(ys_fit)-np.amin(ys_fit))
			for peak in self.peaks:
				x, y, w1, w2 = peak
				if not 0 < w1 < width:
					w1 = width
				if not 0 < w2 < width:
					w2 = width
				if not x_range[0] < x < x_range[1]:
					x = (x_range[0]+x_range[1])/2
				if not 0 < y < amp:
					y = amp
				if self.arg_nr == 4:
					p0.extend((x, y, w1, w2))
					bounds[0].extend((x_range[0], 0, 0, 0))
					bounds[1].extend((x_range[1], amp, width, width))
				else:
					p0.extend((x, y, w1))
					bounds[0].extend((x_range[0], 0, 0))
					bounds[1].extend((x_range[1], amp, width))
			p0.extend([0]*self.polynom)
			bounds[0].extend([-np.inf]*self.polynom)
			bounds[1].extend([+np.inf]*self.polynom)
			popt, pciv = optimize.curve_fit(self.f, self.exp_xs, ys_fit, p0=p0, bounds=bounds)
			ys = self.f(xs, *popt)
		else:
			popt = (0,1,1,1)*len(self.peaks)
			ys = 0*xs
		self._plot_fit[0].set_data(xs, ys)
		opt_param = []
		for i in range(len(self.peaks)):
			ys = self.f(xs, *popt[i*self.arg_nr: (i+1)*self.arg_nr], pol=False)
			tmp_plot = self._ax.plot(xs, ys, color = "blue", alpha=trans)
			self._plots_parts.append(tmp_plot[0])
			opt_param.append( list(popt[i*self.arg_nr: (i+1)*self.arg_nr])+[i]+list(self.peaks[i]) )
		self._plots_parts.append(self._ax.scatter([x[0] for x in opt_param], [x[1] for x in opt_param], color="red"))

		if self.polynom > 0 and self.checkbox_polynom.isChecked():
			self._plots_parts.append(self._ax.plot(xs, np.polyval(popt[-self.polynom:], xs-sum(self.x_range)/2), color="yellow")[0])
		
		opt_param.sort(key=lambda x: x[0])
		self.plot_canvas.draw()
		table = self.table
		table.setRowCount(0)
		table.setColumnCount(10)
		table.setHorizontalHeaderLabels(["Action", "Frequency", "Amplitude", "Width Gauss", "Width Lorentz", "Start_Freq", "Start_Amp", "Start_Gauss", "Start_Lorentz", "Delete"])
		for params in opt_param:
			currRowCount = table.rowCount()
			table.insertRow(currRowCount)
			tmp_button = QPushButton('Assign')
			tmp_button.clicked.connect(lambda a, i=self.i, freq=params[0]: self.main_widget.assign(i, freq=freq))
			table.setCellWidget(currRowCount,0,tmp_button)
			table.setItem(currRowCount, 1, QTableWidgetItem(f'{params[0]:.4f}'))
			table.setItem(currRowCount, 2, QTableWidgetItem(f'{params[1]:.4f}'))
			table.setItem(currRowCount, 3, QTableWidgetItem(f'{params[2]:.4f}'))
			table.setItem(currRowCount, 4, QTableWidgetItem(f'{params[3]:.4f}'))
			table.setItem(currRowCount, 5, QTableWidgetItem(f'{params[5]:.4f}'))
			table.setItem(currRowCount, 6, QTableWidgetItem(f'{params[6]:.4f}'))
			table.setItem(currRowCount, 7, QTableWidgetItem(f'{params[7]:.4f}'))
			table.setItem(currRowCount, 8, QTableWidgetItem(f'{params[8]:.4f}'))
			tmp_button2 = QPushButton('Delete')
			tmp_button2.clicked.connect(lambda a, ind=params[4]: self.del_peak(i=ind))
			table.setCellWidget(currRowCount,9,tmp_button2)
		table.resizeColumnsToContents()

	@synchronized(locks["windows"])
	def closeEvent(self, event):
		self.config["blendedlineswindow_lineshape"] = self.input_func.currentText()
		self.config["blendedlineswindow_derivative"] = self.input_derivative.text()
		self.config["blendedlineswindow_transparency"] = self.input_transp.value()
		self.config["blendedlineswindow_maxfwhm"] = self.input_width.text()
		self.config["blendedlineswindow_polynom"] = self.input_polynom.value()
		self.config["blendedlineswindow_showbaseline"] = self.checkbox_polynom.isChecked()
		self.fig.canvas.mpl_disconnect(self.cid)
		super().closeEvent(event)
		
##
## Propose Series to start with
##
class PathfinderWindow(EQWidget):
	def __init__(self, parent=None):
		super(PathfinderWindow, self).__init__(parent)
		self.setWindowTitle("Pathfinder")
		self.resize(1000, 450)

		self.main_window = mw
		self.main_widget = mw.main_widget
		self.config = self.main_widget.config
		layout = QVBoxLayout()

		introductionLabel = QLabel("Welcome to Pathfinder,\nhere we help you to get a foot in the door.")
		introductionLabel.setWordWrap(True)

		rangeLabel = QLabel("Frequency Range:")
		rangeStartLabel = QLabel("Start Frequency: ")
		self.rangeStartInput = QLineEdit()
		rangeStopLabel = QLabel("Stop Frequency: ")
		self.rangeStopInput = QLineEdit()

		self.startJourneyButton = QPushButton("Go!")
		self.startJourneyButton.clicked.connect(lambda a: self.show_largest())

		numberOfResultsLabel = QLabel("Number of Results: ")
		self.numberOfResults = QSpinBox()
		self.numberOfResults.setMaximum(100)
		self.numberOfResults.setValue(10)

		additionalConditionLabel = QLabel("Additional Condition: ")
		self.addConditionInput = QLineEdit()

		self.messageLabel = QLabel()
		self.messageLabel.setWordWrap(True)
		self.messageLabel.setHidden(True)
		self.messageLabel.setHidden(True)

		self.outputTable = QTableWidget()
		self.outputTable.setEditTriggers(QAbstractItemView.NoEditTriggers)
		self.outputTable.setHidden(True)
		self.outputTable.setMinimumHeight(200)
		self.outputTable.move(0,0)

		self.transitions_type_cbs = []
		transitionsLabel = QLabel("Allowed Transitions")
		tmp = QCheckBox("a-type")
		self.transitions_type_cbs.append(tmp)
		tmp = QCheckBox("b-type")
		self.transitions_type_cbs.append(tmp)
		tmp = QCheckBox("c-type")
		self.transitions_type_cbs.append(tmp)

		vertLayout = QHBoxLayout()
		leftLayout = QGridLayout()
		rightLayout = QVBoxLayout()
		layout.addWidget(introductionLabel)

		rightLayout.addWidget(transitionsLabel)
		rightLayout.addWidget(self.transitions_type_cbs[0])
		rightLayout.addWidget(self.transitions_type_cbs[1])
		rightLayout.addWidget(self.transitions_type_cbs[2])
		rightLayout.addStretch(1)


		i=0
		leftLayout.addWidget(rangeLabel, i, 0)
		i+=1
		leftLayout.addWidget(rangeStartLabel, i, 0)
		leftLayout.addWidget(self.rangeStartInput, i, 1)
		i+=1
		leftLayout.addWidget(rangeStopLabel, i, 0, 1, 1)
		leftLayout.addWidget(self.rangeStopInput, i, 1, 1, 1)
		i+=1
		leftLayout.addWidget(numberOfResultsLabel, i, 0, 1, 1)
		leftLayout.addWidget(self.numberOfResults, i, 1, 1, 1)
		i+=1
		leftLayout.addWidget(additionalConditionLabel, i, 0, 1, 1)
		leftLayout.addWidget(self.addConditionInput, i, 1, 1, 1)
		i+=1
		leftLayout.addWidget(self.startJourneyButton, i, 1, 1, 1)
		i+=1

		vertLayout.addLayout(leftLayout)
		vertLayout.addLayout(rightLayout)
		vertLayout.addStretch(0)
		layout.addLayout(vertLayout)
		layout.addWidget(self.messageLabel)
		layout.addWidget(self.outputTable)
		layout.addStretch(1)
		self.setLayout(layout)
		
		self.rangeStartInput.setText(self.config["pathfinderwindow_start"])
		self.rangeStopInput.setText(self.config["pathfinderwindow_stop"])
		self.numberOfResults.setValue(self.config["pathfinderwindow_results"])
		self.addConditionInput.setText(self.config["pathfinderwindow_condition"])
		[x.setChecked(True) for x,y in zip(self.transitions_type_cbs, json.loads(self.config["pathfinderwindow_types"])) if y]

	@synchronized(locks["cat_df"])
	@synchronized(locks["ass_df"])
	def show_largest(self):
		condition = []
		nor = int(self.numberOfResults.value())
		transition_cbs = [x.isChecked() for x in self.transitions_type_cbs]
		tmp_min = self.rangeStartInput.text()
		tmp_max = self.rangeStopInput.text()
		addCondition = self.addConditionInput.text()
		if addCondition.strip()!="":
			condition.append(addCondition)
		try:
			tmp_min = float(tmp_min)
		except ValueError:
			tmp_min = None
		try:
			tmp_max = float(tmp_max)
		except ValueError:
			tmp_max = None
		if tmp_min != None:
			condition.append(f"{tmp_min} <= x")
		if tmp_max != None:
			condition.append(f"x <= {tmp_max}")

		tmp_condition = []
		if transition_cbs[0] == True:
			tmp_condition.append(f"abs(qnu2-qnl2) == 0 and abs(qnu3-qnl3) == 1")
		if transition_cbs[1] == True:
			tmp_condition.append(f"abs(qnu2-qnl2) == 1 and abs(qnu3-qnl3) == 1")
		if transition_cbs[2] == True:
			tmp_condition.append(f"abs(qnu2-qnl2) == 1 and abs(qnu3-qnl3) == 0")
		if tmp_condition != []:
			condition.append(" or ".join([f"({x})" for x  in tmp_condition]))
		tmp_df = self.main_widget._lin_data_df.copy()
		querystring = " and ".join([f"({x})" for x  in condition])
		if querystring != "":
			tmp_df.query(querystring, inplace=True)
		if self.config["flag_debug"]:
			self.main_widget.notification(f"Querystring is {querystring}.")
			self.main_widget.notification(f"Found {len(tmp_df)} results matching the query.")
		all_qns_labels = ["qnu1", "qnu2", "qnu3", "qnu4", "qnu5", "qnu6", "qnl1", "qnl2", "qnl3", "qnl4", "qnl5", "qnl6"]
		visible = [not x.isHidden() for x in self.main_widget.QNcbs]
		qns_labels = [x for i,x in enumerate(all_qns_labels) if visible[i%len(visible)]]
				
		tmp_catass_df = self.main_widget._catass_df.copy()
		tmp_catass_df["DROP"] = True
		tmp_catass_df.drop(columns=["x", "y"]+list(set(all_qns_labels)-set(qns_labels)), inplace=True)
		tmp_df = pd.merge(tmp_df, tmp_catass_df, how="outer", on=qns_labels)
		tmp_df = tmp_df[tmp_df.DROP != True]

		largest = tmp_df.nlargest(nor, "y")

		if tmp_min != None and tmp_max != None:
			range = f" in the range from {tmp_min:g} to {tmp_max:g}"
		else:
			range = " without a specified range"
		message = f"We found some some very promising starting points{range}. Here are the {nor} best starting points, according to the intensities given in the .cat file. Already assigned lines are excluded."
		self.messageLabel.setText(message)
		self.messageLabel.setHidden(False)

		table = self.outputTable
		headers = ["Start", "Intensity", "Frequency", "qnu1", "qnu2", "qnu3", "qnu4", "qnu5", "qnu6", "qnl1", "qnl2", "qnl3", "qnl4", "qnl5", "qnl6"]
		table.setRowCount(0)
		table.setColumnCount(len(headers))
		table.setHorizontalHeaderLabels(headers)
		new_headers = []
		for index, row in largest.iterrows():
			currRowCount = table.rowCount()
			table.insertRow(currRowCount)
			table.setItem(currRowCount, 1, QTableWidgetItem(f'{row["y"]:.4f}'))
			table.setItem(currRowCount, 2, QTableWidgetItem(f'{row["x"]:.4f}'))
			tmp_button = QPushButton('Start')
			tmp_button.clicked.connect(lambda a,crow=row: self.startHere(crow))
			table.setCellWidget(currRowCount,0,tmp_button)
			i = 3
			for column in ["qnu1", "qnu2", "qnu3", "qnu4", "qnu5", "qnu6"]:
				if row[column] == -1:
					break
				table.setItem(currRowCount, i, QTableWidgetItem(f'{row[column]:g}'))
				new_headers.append(column)
				i+=1
			for column in ["qnl1", "qnl2", "qnl3", "qnl4", "qnl5", "qnl6"]:
				if row[column] == -1:
					break
				table.setItem(currRowCount, i, QTableWidgetItem(f'{row[column]:g}'))
				new_headers.append(column)
				i+=1
		new_headers = [x for x in headers if x in list(set(new_headers))]
		new_headers = new_headers if len(new_headers) != 0 else ["qnu1", "qnu2", "qnu3", "qnl1", "qnl2", "qnl3"]
		headers = ["Start", "Intensity", "Frequency"]+new_headers
		table.setColumnCount(len(headers))
		table.setHorizontalHeaderLabels(headers)
		self.outputTable.resizeColumnsToContents()
		self.outputTable.setHidden(False)

	def startHere(self, row):
		QNUs = [row["qnu1"], row["qnu2"], row["qnu3"], row["qnu4"], row["qnu5"], row["qnu6"]]
		QNLs = [row["qnl1"], row["qnl2"], row["qnl3"], row["qnl4"], row["qnl5"], row["qnl6"]]
		DIFF = [QNUs[i]-QNLs[i] for i in range(len(QNUs))]
		for i in range(len(QNUs)):
			if QNUs[i] == -1 or QNLs[i]  == -1:
				self.main_widget.alter_QNs(absolute=i)
				break
			else:
				self.main_widget.QNUs[i].setValue(QNUs[i])
				self.main_widget.QNLs[i].setValue(QNLs[i])
				if DIFF[i] != 0:
					self.main_widget.QNcbs[i].setChecked(True)
		self.main_widget.tabs.setCurrentIndex(1)
		self.main_widget.update_plots()
		self.main_widget.activateWindow()

	@synchronized(locks["windows"])
	def closeEvent(self, event):
		self.config["pathfinderwindow_start"]		=	self.rangeStartInput.text()
		self.config["pathfinderwindow_stop"]		=	self.rangeStopInput.text()
		self.config["pathfinderwindow_results"]		=	self.numberOfResults.value()
		self.config["pathfinderwindow_condition"]	=	self.addConditionInput.text()
		self.config["pathfinderwindow_types"]		=	json.dumps([x.isChecked() for x in self.transitions_type_cbs])
		super().closeEvent(event)

##
## Allow for setting up and starting Pipe
##
class PipeWindow(EQWidget):
	def __init__(self, parent=None):
		super(PipeWindow, self).__init__(parent)
		self.setWindowTitle("Pipe")

		self.main_window = mw
		self.main_widget = mw.main_widget
		self.config = self.main_widget.config

		self.welcomeLabel = QLabel("Here you can set up a pipe to run fitting and prediction programs to update your files.")
		self.welcomeLabel.setWordWrap(True)

		self.commandInput = QPlainTextEdit(self.config["pipe_command"])
		self.commandInput.textChanged.connect(lambda text="": self.changedLabel.setText("Unsaved Changes"))

		self.changedLabel = QLabel("")

		self.rereadCheckbox = QCheckBox("Reread Files after command finished")
		if self.config["pipe_reread"] == True:
			self.rereadCheckbox.setChecked(True)
		self.rereadCheckbox.stateChanged.connect(lambda x:self.update_reread())

		self.saveButton = QPushButton("Save")
		self.saveButton.clicked.connect(self.update_pipe_command)

		self.runButton = QPushButton("Run")
		self.runButton.clicked.connect(lambda x: self.main_window.run_pipe())


		layout = QVBoxLayout()
		layout.addWidget(self.welcomeLabel)
		layout.addWidget(self.commandInput)
		layout.addWidget(self.changedLabel)
		layout.addWidget(self.rereadCheckbox)
		layout.addWidget(self.saveButton)
		layout.addWidget(self.runButton)
		self.setLayout(layout)

	def update_pipe_command(self):
		command = self.commandInput.toPlainText()
		if command.strip() == "":
			command = None
		self.config["pipe_command"] = command
		self.changedLabel.setText("")

	def update_reread(self):
		self.config["pipe_reread"] = self.rereadCheckbox.isChecked()

##
## Controls for Fitting Series
##
class SeriesFitWindow(EQWidget):
	def __init__(self, parent=None):
		super(SeriesFitWindow, self).__init__(parent)
		self.setWindowTitle("Series Fit")

		self.main_window = mw
		self.main_widget = mw.main_widget
		self.config = self.main_widget.config
		
		self.pred_function = None
		self.pred_params = {}

		self.catalogueTableDeleteButton = QPushButton("X")
		self.catalogueTableDeleteButton.setDefault(True)
		self.catalogueTableDeleteButton.clicked.connect(lambda x: self.catalogueTableDelete())
		self.catalogueTableDeleteButton.setFixedWidth(20)

		self.catalogueTableResizeButton = QPushButton("Tight")
		self.catalogueTableResizeButton.setDefault(True)
		self.catalogueTableResizeButton.clicked.connect(lambda x: self.catalogueTable.resizeColumnsToContents())
		self.catalogueTableResizeButton.setFixedWidth(50)

		self.catalogueTable = QTableView()
		self.catalogueModel = CatalogueTableModel(self.main_widget._seriesfit, self.main_widget._seriesfit_columns)
		self.catalogueTable.setModel(self.catalogueModel)
		
		
		tab2layout = QGridLayout()
		self.QNU1 = QSpinBox()
		self.QNU2 = QSpinBox()
		self.QNU3 = QSpinBox()
		self.QNU4 = QSpinBox()
		self.QNU5 = QSpinBox()
		self.QNU6 = QSpinBox()

		self.QNUs = [self.QNU1, self.QNU2, self.QNU3, self.QNU4, self.QNU5, self.QNU6]

		self.QNL1 = QSpinBox()
		self.QNL2 = QSpinBox()
		self.QNL3 = QSpinBox()
		self.QNL4 = QSpinBox()
		self.QNL5 = QSpinBox()
		self.QNL6 = QSpinBox()

		self.QNLs = [self.QNL1, self.QNL2, self.QNL3, self.QNL4, self.QNL5, self.QNL6]

		self.QNcb1 = QCheckBox("Incr")
		self.QNcb2 = QCheckBox("Incr")
		self.QNcb3 = QCheckBox("Incr")
		self.QNcb4 = QCheckBox("Incr")
		self.QNcb5 = QCheckBox("Incr")
		self.QNcb6 = QCheckBox("Incr")

		self.QNcbs = [self.QNcb1, self.QNcb2, self.QNcb3, self.QNcb4, self.QNcb5, self.QNcb6]

		self.IncQN = QPushButton("Incr")
		self.DecQN = QPushButton("Decr")
		self.IncQN.clicked.connect(lambda x: self.QNs_incdec(+1))
		self.DecQN.clicked.connect(lambda x: self.QNs_incdec(-1))

		self.IncNumberQNs = QPushButton("+")
		self.DecNumberQNs = QPushButton("-")
		self.IncNumberQNs.clicked.connect(lambda x: self.alter_QNs(+1))
		self.DecNumberQNs.clicked.connect(lambda x: self.alter_QNs(-1))

		self.labels = (QLabel("QN 1"), QLabel("QN 2"), QLabel("QN 3"), QLabel("QN 4"), QLabel("QN 5"), QLabel("QN 6"))
		for i in range(len(self.labels)):
			tab2layout.addWidget(self.labels[i], 1, i+1)

		for i in range(len(self.QNUs)):
			tab2layout.addWidget(self.QNUs[i], 2, i+1)
			self.QNUs[i].setMaximumWidth(50)
			self.QNUs[i].setMaximum(10000)

		for i in range(len(self.QNLs)):
			tab2layout.addWidget(self.QNLs[i], 3, i+1)
			self.QNLs[i].setMaximumWidth(50)
			self.QNLs[i].setMaximum(10000)

		for i in range(len(self.QNcbs)):
			tab2layout.addWidget(self.QNcbs[i], 4, i+1)
			self.QNcbs[i].setMaximumWidth(50)

		tab2layout.addWidget(self.IncQN, 2, 7, 1, 2)
		tab2layout.addWidget(self.DecQN, 3, 7, 1, 2)

		tab2layout.addWidget(self.IncNumberQNs, 1, 7)
		tab2layout.addWidget(self.DecNumberQNs, 1, 8)

		tab2layout.setColumnStretch(10, 1)
		self.fitInput = QPlainTextEdit()
		self.fitInput.setMinimumHeight(50)
		self.fitButton = QPushButton("Fit")
		self.fitButton.clicked.connect(self.fit)
		self.fitButton.setMaximumWidth(100)
		self.textArea = QPlainTextEdit()
		self.textArea.setReadOnly(True)
		self.textArea.setMinimumHeight(50)

		layout = QVBoxLayout()
		buttonsBox = QHBoxLayout()
		buttonsBox.addWidget(self.catalogueTableDeleteButton)
		buttonsBox.addStretch(2)
		buttonsBox.addWidget(self.catalogueTableResizeButton)
		layout.addLayout(buttonsBox)
		layout.addWidget(self.catalogueTable)
		layout.addLayout(tab2layout)
		layout.addWidget(self.fitInput)
		layout.addWidget(self.fitButton)
		layout.addWidget(self.textArea)
		self.setLayout(layout)
		
		self.alter_QNs(0)
		if self.config["seriesfitwindow_transition"] != None:
			[inp.setValue(val) for val, inp in zip(self.config["series_transition"][0], self.QNUs)]
			[inp.setValue(val) for val, inp in zip(self.config["series_transition"][1], self.QNLs)]
			[cb.setChecked(val) for val, cb in zip(self.config["series_transition"][2], self.QNcbs)]
		self.fitInput.document().setPlainText(self.config["seriesfitwindow_function"])
	
	def alter_QNs(self, index=None, absolute=None):
		if index != None:
			self.config["seriesfitwindow_qns"] += index
		elif absolute != None:
			self.config["seriesfitwindow_qns"] = absolute
		if self.config["seriesfitwindow_qns"] > 6:
			self.config["seriesfitwindow_qns"] = 6
		elif self.config["seriesfitwindow_qns"] < 1:
			self.config["seriesfitwindow_qns"] = 1
		for i in range(self.config["seriesfitwindow_qns"]):
			self.QNUs[i].setHidden(False)
			self.QNLs[i].setHidden(False)
			self.QNcbs[i].setHidden(False)
			self.labels[i].setHidden(False)

		for i in range(self.config["seriesfitwindow_qns"], len(self.labels)):
			self.QNUs[i].setHidden(True)
			self.QNLs[i].setHidden(True)
			self.QNcbs[i].setHidden(True)
			self.labels[i].setHidden(True)

	@synchronized(locks["ser_df"])
	def catalogueTableDelete(self):
		selected = [index.row() for index in self.catalogueTable.selectionModel().selectedRows()]
		for index in sorted(selected, reverse = True):
			self.main_widget._seriesfit.drop(index, inplace =True)
			self.main_widget._seriesfit.reset_index(inplace = True, drop = True)
		self.catalogueModel.update()
		
	def fit(self):
		qnu = [float(x.text()) if x.text()!= "" else 0 for x in self.QNUs]
		qnl = [float(x.text()) if x.text()!= "" else 0 for x in self.QNLs]
		incr = [1 if x.isChecked() else 0 for x in self.QNcbs]
		diff = [1,1,1,1,1,1]
		visible = [not x.isHidden() for x in self.QNcbs]
		
		labels_upper = ["qnu1", "qnu2", "qnu3", "qnu4", "qnu5", "qnu6"]
		labels_lower = ["qnl1", "qnl2", "qnl3", "qnl4", "qnl5", "qnl6"]
		
		df = self.main_widget._seriesfit.copy()
		
		condition_incr = []
		condition_noincr = []
		for i in range(len(incr)):
			if visible[i]:
				if incr[i]:
					condition_incr.append(f"{labels_upper[i]} - {qnu[i]}")
					condition_incr.append(f"{labels_lower[i]} - {qnl[i]}")
				else:
					condition_noincr.append(f"{labels_upper[i]} == {qnu[i]}")
					condition_noincr.append(f"{labels_lower[i]} == {qnl[i]}")
		
		conditions = []
		if len(condition_incr) != 0:
			conditions.append(" == ".join(condition_incr))
		conditions.extend(condition_noincr)
		query_string = " & ".join(conditions)
		
		tmp_df = df.query(query_string).copy()
		tmp_df.sort_values("x", inplace = True)
		
		eval_str = self.fitInput.toPlainText()
		if eval_str == "":
			eval_str = "0"
		
		self.fit_code = compile(eval_str, "<string>", "eval")
		self.fit_params = set(self.fit_code.co_names)-set(["qnu1", "qnu2", "qnu3", "qnu4", "qnu5", "qnu6","qnl1", "qnl2", "qnl3", "qnl4", "qnl5", "qnl6"])
		
		if len(tmp_df) < len(self.fit_params):
			self.writeLog("Too few assignments for too many parameters.\n\n")
			return
		elif len(self.fit_params) == 0:
			self.writeLog("No parameters in your expression.\n\n")
			return
			
		
		p0 = [1]*len(self.fit_params)
		popt, pciv = optimize.curve_fit(self.function, tmp_df[labels_upper+labels_lower].values.transpose(), tmp_df["x"].to_numpy(), p0=p0)
		
		qnu = np.array(qnu)
		qnl = np.array(qnl)
		incr = np.array(incr)
		qnss = [np.concatenate((qnu+incr*i, qnl+incr*i)) for i in range(100)]
		

		freqs = [self.function(qns, *popt) for qns in qnss]
		self.main_widget.tabs.setCurrentIndex(2)
		self.main_widget.load_frequencies_list(freqs)
		
		tmp = [f"{name} : {value:.4f}" for name, value in zip(self.fit_params, popt)]
		tmp = "\n".join(tmp)
		self.writeLog(f"Succeeded, parameters were determined as \n{tmp}")
		
	def writeLog(self, text):
		time_str = time.strftime("%H:%M", time.localtime())
		self.textArea.appendPlainText(f"{time_str}: {text}")
		sb = self.textArea.verticalScrollBar()
		sb.setValue(sb.maximum())
	
	def function(self, qns, *params):
		qnu1, qnu2, qnu3, qnu4, qnu5, qnu6, qnl1, qnl2, qnl3, qnl4, qnl5, qnl6 = qns
		tmp_dict = {
			"qnu1" : qnu1,
			"qnu2" : qnu2,
			"qnu3" : qnu3,
			"qnu4" : qnu4,
			"qnu5" : qnu5,
			"qnu6" : qnu6,
			"qnl1" : qnl1,
			"qnl2" : qnl2,
			"qnl3" : qnl3,
			"qnl4" : qnl4,
			"qnl5" : qnl5,
			"qnl6" : qnl6,
		}
		for name, param in zip(self.fit_params, params):
			tmp_dict[name] = param
			
		return(eval(self.fit_code, tmp_dict))
	
	@synchronized(locks["windows"])
	def closeEvent(self, event):
		self.config["seriesfitwindow_transition"]		= json.dumps([	[x.value() for x in self.QNUs],
																[x.value() for x in self.QNLs],
																[x.isChecked() for x in self.QNcbs],
															])
		self.config["seriesfitwindow_function"] 		= self.fitInput.toPlainText()
		super().closeEvent(event)

##
## Credits Window
##
class CreditsWindow(EQWidget):
	def __init__(self, parent=None):
		super(CreditsWindow, self).__init__(parent)
		self.setWindowTitle("Credits")

		self.main_window = mw
		self.main_widget = mw.main_widget
		self.config = self.main_widget.config

		layout = QVBoxLayout()

		global credits_string
		self.textWidget = QPlainTextEdit(credits_string)
		self.textWidget.setReadOnly(True)
		self.textWidget.setMinimumHeight(300)
		self.textWidget.setMinimumWidth(300)
		layout.addWidget(self.textWidget)

		self.setLayout(layout)
		

##
## SignalClass for Notifications from non main-thread and Error Class for Custom Error
##
class SignalClass(QObject):
	writeLog = pyqtSignal(str)
	updateWindows = pyqtSignal()
	
	def __init__(self):
		super(SignalClass, self).__init__()

class CustomError(Exception):
	pass

##
## Pandas Table Model
##
class CatalogueTableModel(QAbstractTableModel):

	def __init__(self, data, headers=[], invisible_trailing=0):
		super(CatalogueTableModel, self).__init__()
		self._data = data
		self._headers = headers
		self.invisible_trailing = invisible_trailing
		
	def data(self, index, role):
		if role == Qt.DisplayRole or role == Qt.EditRole:
			value = self._data.iloc[index.row(), index.column()]
			if value == -1:
				return("")
			elif int(value) == value:
				return(f"{value:.0f}")
			else:
				return (f"{value:.4f}")

	def rowCount(self, index):
		return self._data.shape[0]

	def columnCount(self, index):
		#-2 for two hidden columns freq_predicted and intens_predicted
		return self._data.shape[1]-self.invisible_trailing
	
	def headerData(self, section, orientation, role):
		# section is the index of the column/row.
		if role == Qt.DisplayRole:
			if orientation == Qt.Horizontal:
				if section in range (len(self._headers)):
					return(str(self._headers[section]))
				else:
					return str(self._data.columns[section])

			if orientation == Qt.Vertical:
				return str(self._data.index[section])
	
	def flags(self, index):
		return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable
	
	def update(self):
		self.layoutChanged.emit()
	
	def setData(self, index, value, role):
		try:
			value = np.float64(value)
		except ValueError:
			value = -1
		if not index.isValid():
			return False
		if role != Qt.EditRole:
			return False
		row = index.row()
		if row < 0 or row >= len(self._data.values):
			return False
		column = index.column()
		if column < 0 or column >= self._data.columns.size:
			return False
		self._data.values[row][column] = value
		self.dataChanged.emit(index, index)
		return True 

##
## Global Function
##
def match_color(color):
	match = re.search(r'^#(?:[0-9a-fA-F]{6})$|^#(?:[0-9a-fA-F]{8})$', color)
	if match:
		return(color)
	else:
		return(False)

def lineshape(shape, derivative, *args):
	if shape == "Gauss":
		x, x_0, amp, width = args
		width = width/(2*np.sqrt(2*np.log(2)))
		ys = np.exp(-(x-x_0)**2/(2*width**2))/(width*(2*np.pi)**0.5)
	elif shape == "Lorentz":
		x, x_0, amp, width = args
		width = width/2
		ys = 1/(np.pi*width*(1+((x-x_0)/width)**2))
	elif shape == "Voigt":
		x, x_0, amp, gauss, lorentz = args
		gauss = gauss/(2*np.sqrt(2*np.log(2))) / np.sqrt(2)
		lorentz = lorentz/2 / np.sqrt(2)
		ys = special.voigt_profile(x-x_0, gauss, lorentz)

	for i in range(0, int(derivative)):
		ys = np.gradient(ys)
	if derivative%2 == 0:
		ys = -ys
	ymax = max(ys)
	if not np.isfinite(ymax):
		ymax = 1
	ys = amp*ys/ymax
	return(ys)


# Generate with following command, quotes is str of quotes.txt
# json.dumps(quotes.split("\n\n"))
quotes_str = '["Deine Leute alles Platzpatronen,\\nMeine Leute halten stand wie ein Bataillon\\n-Anis M. Y. Ferchichi", "Was sie rappen ist kein Image, das ist Science Fiction\\n-Anis M. Y. Ferchichi", "Erstens kommt es anders und zweitens als man denkt\\n-Ekrem Bora", "In der Theorie willst du was ver\\u00e4ndern, aber in der Praxis geht alles schief.\\nIn der Theorie wollen wir zu viel, geben aber viel zu wenig, das ist kein fairer Deal.\\n-Anis M. Y. Ferchichi & Jochen Burchard", "Und bist du unten, dr\\u00fccken sie dich noch ein St\\u00fcck tiefer\\nNoch ein St\\u00fcck tiefer\\nNoch ein St\\u00fcck tiefer\\nNoch ein St\\u00fcck tiefer\\n-Anis M. Y. Ferchichi", "Da ting goes skrrrrahh, pap, pap, ka-ka-ka\\nSkidiki-pap-pap, and a pu-pu-pudrrrr-boom\\nSkya, du-du-ku-ku-dun-dun\\nPoom, poom, you dun know\\n-Michael Dapaah", "Einstein, my mind, this MC move white squares with my relatives.\\n-Shawn C. Carter", "Man kann es weder sehen noch ber\\u00fchren, lediglich sp\\u00fcren.\\n-Anis M. Y. Ferchichi", "Dank mir ist das Niveau im Keller, so wie Weinf\\u00e4sser\\n-Anis M. Y. Ferchichi", "Meine Punchlines sind w\\u00e4hlerisch, aber dich zers\\u00e4ge ich\\n-Vladislav Balovatsky", "Bratan l\\u00e4uft, Para und Erfolg\\nGestern war\'n die Platten grau, heute sind sie Gold\\n-Vladislav Balovatsky", "Das Leben ist \'ne Schachtel Pralinen\\nUnd alle kriegen das, was sie verdienen\\nEgal, wie schnell du rennst, vor Gott kannst du nicht fliehen\\nInshallah deine S\\u00fcnden werden verziehen\\n-Vladislav Balovatsky", "F\\u00fcr meine Tochter, der ich sp\\u00e4ter mal erz\\u00e4hlen kann\\nIhr Vater war ein Ehrenmann\\n-Anis M. Y. Ferchichi", "Nachts, da schlafe ich\\n-Unbekannter deutscher Philosoph", "It wasn\'t me\\n-Orville Richard Burrell", "Und du bist, alles f\\u00fcr mich, alles das, was mir Angst macht\\nMeine silberne Kugel, mein Kryptonit und mein Anthrax\\n-Friedrich Kautz", "Jack-Daniel\'s-Eis, Hennessy und Johnnie Walker\\nIch komm\' voll besoffen in \'nem Gucci-Jogger\\n-Vladislav Balovatsky", "There is no such thing as too much duct tape\\n-Oli", "We are going so fast\\nBut time so slow\\n-Theory of Relativity", "D\\u00f6p, d\\u00f6p, d\\u00f6p, d\\u00f6p, d\\u00f6d\\u00f6d\\u00f6d\\u00f6d\\u00f6p\\nD\\u00f6, d\\u00f6d\\u00f6d\\u00f6p, d\\u00f6p, d\\u00f6d\\u00f6d\\u00f6d\\u00f6d\\u00f6p, d\\u00f6p\\n-Hans Peter Geerdes", "Skibadee, skibadanger\\nI am the rearranger\\n-Hans Peter Geerdes", "If we could only slow the time\\nWe would have forever every night\\n-Don Pepijn Schipper", "Meine Stra\\u00dfenpoesie l\\u00f6st die Chaostheorie\\n-Mousa Amouei", "Die Parabel sie steigt, und zwar exponentiell\\n-Mohamed El Moussaoui", "Chuba chuba chuba chuba chuba chuba chubby.\\nI don\'t have any lines to go right here, so chuby Teletubby\\n-Marshall Bruce Mathers III", "Two things are infinite: The universe and human stupidity;\\nand I\\u2018m not sure about the universe\\n-Albert E", "Physics is like sex: sure, it may give some practical results, but that\'s not why we do it.\\n-Richard P. Feynman", "I do not think you can name many great inventions that have been made by married men.\\n-Nikola Tesla", "Those who are not shocked when they first come across quantum theory cannot possibly have understood it.\\n-Niels Bohr", "We\\u2019re not free in what we do, because we\\u2019re not free in what we want.\\n-Jonas Kahnwald", "What we know is a drop. What we don\\u2019t know is an ocean.\\n-Isaac Newton", "Oh fate, you mysterious bitch\\nWhy you doing this to me\\n-Theodore \\"T-Bag\\" Bagwell", "We are captives of our own identities, living in prisons of our own creation.\\n-Theodore \\"T-Bag\\" Bagwell", "Preparation can only take you so far. \\nAfter that, you have to take a few leaps of faith.\\n-Michael Scofield", "Das ist der Sound f\\u00fcr die echten M\\u00e4nner, die das hier h\\u00f6ren, wenn sie Pressluft h\\u00e4mmern\\n-Tarek Ebene, Nico Seyfrid & Maxim Dr\\u00fcner"]'
def except_hook(cls, exception, traceback):
	sys.__excepthook__(cls, exception, traceback)
	with open(f"{app_tag}.err", "a+") as file:
		time_str = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
		file.write(f"{time_str}: \n{exception}\n{''.join(tb.format_tb(traceback))}\n\n")
	try:
		if not mw.main_widget._catalogue.empty:
			mw.save_lines_cat(f"{app_tag}.cat", force_new=True)
	except Exception as E:
		pass


if __name__ == '__main__':
	sys.excepthook = except_hook
	app = QApplication(sys.argv)

	mw = MainWindow()
	mw.show()

	sys.exit(app.exec_())
