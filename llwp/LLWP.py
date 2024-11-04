#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Luis Bonah
# Description: Loomis-Wood Plot Software for Assigning experimental Spectra to Quantum Numbers

CREDITSSTRING = """Made by Luis Bonah

As this programs GUI is based on PyQt6, which is GNU GPL v3 licensed, this program is also licensed under GNU GPL v3 (See the bottom paragraph).

pandas, matplotlib, scipy and numpy were used for this program, speeding up the development process massively.

Copyright (C) 2024

	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program. If not, see <https://www.gnu.org/licenses/>.
"""

##
## Global Constants and Imports
##

APP_TAG = "LLWP"

import os
import sys
import re
import io
import csv
import time
import wrapt
import pickle
import json
import threading
import configparser
import traceback as tb
import numpy as np
import pandas as pd
import subprocess
import webbrowser
import pyckett

from functools import lru_cache
from scipy import optimize, special, signal

from PyQt6.QtCore import *
from PyQt6.QtWidgets import *
from PyQt6.QtGui import *

import matplotlib
from matplotlib.backends.backend_qtagg import FigureCanvas, NavigationToolbar2QT

import warnings

try:
	warnings.simplefilter('ignore', np.RankWarning)
except AttributeError:
	warnings.simplefilter('ignore', np.exceptions.RankWarning)

# Can be changed right before starting program:
# On Unix systems:
# export LLWP_QNS=10
# On Windows systems:
# set LLWP_QNS=10
N_QNS = int(os.environ.get('LLWP_QNS', 10))
N_QNS = max(6, N_QNS)

QLocale.setDefault(QLocale('en_EN'))
matplotlib.rcParams['axes.formatter.useoffset'] = False

def status_d(func):
	def _wrapper(*args, **kwargs):
		mw = mainwindow
		tmp = mw.status_counter.increase()
		if tmp:
			mw.working.emit()

		try:
			return(func(*args, **kwargs))
		except Exception as E:
			raise
		finally:
			tmp = mw.status_counter.decrease()
			if not tmp:
				mw.waiting.emit()
	return(_wrapper)

def lock_d(lock):
	@wrapt.decorator
	def _wrapper(wrapped, instance, args, kwargs):
		with lock:
			return wrapped(*args, **kwargs)
	return _wrapper

class Color(str):
	@staticmethod
	def trgb_to_rgbt(color):
		if len(color) == 9:
			color = f"#{color[3:]}{color[1:3]}"
		return(color)

	@staticmethod
	def rgbt_to_trgb(color):
		if len(color) == 9:
			color = f"#{color[-2:]}{color[1:-2]}"
		return(color)

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

class CustomError(Exception):
	pass

class EarlyReturnError(Exception):
	pass

class GUIAbortedError(Exception):
	pass

class AtomicCounter():
    def __init__(self, initial=0, callback=None):
        self.value = initial
        self._lock = threading.Lock()
        
        if callback is not None:
            self.callback = callback

    def increase(self, num=1):
        with self._lock:
            self.value += num
            return self.value

    def decrease(self, num=1):
        return(self.increase(-num))

    def get_value(self):
        with self._lock:
            return(self.value)

    def __enter__(self, *args):
        self.increase()

    def __exit__(self, type, value, traceback):
        if self.decrease() == 0:
            self.callback()

    def callback(self):
        pass

class DynamicDecorator():
	def __init__(self): 
		self.init_func = lambda *args, **kwargs: None
		self.exit_func = lambda *args, **kwargs: None
		self.callback_kwargs = {}

	def d(self, func):
		def _wrapper(*args, **kwargs):
			self.init_func(**self.callback_kwargs)
			try:
				return(func(*args, **kwargs))
			except Exception as E:
				raise E
			finally:
				self.exit_func(**self.callback_kwargs)
		return(_wrapper)	

class Config(dict):
	initial_values = {
		'plot_dpi': (100, int),
		'plot_width': (20, float),
		'plot_offset': (0, float),
		'plot_offsetisrelative': (True, bool),
		'plot_rows': (5, int),
		'plot_cols': (1, int),
		'plot_widthexpression': ('', str),
		'plot_offsetexpression': ('', str),

		'plot_yscale': ("Per Plot", str),
		'plot_expcat_factor': (1, float),
		'plot_expcat_exponent': (10, int),
		'plot_yscale_min': (-100, float),
		'plot_yscale_max': (300, float),
		'plot_ymargin': (0.1, float),
		'plot_xticks': (3, int),
		'plot_xtickformat': ('offset', str),

		'plot_gridspeckwargs': ({"hspace": 0, "wspace": 0}, dict),
		'plot_fontdict': ({"size":10}, dict),
		'plot_bins': (4000, int),
		'plot_skipbinning': (1000, int),
		'plot_hovercutoff': (20, float),
		'plot_annotationkwargs': ({"x": 1, "y": 1, "horizontalalignment": "right", "verticalalignment": "top"}, dict),
		'plot_annotationfstring': ('{x}', str),

		'convolution_kernel': ([], list),
		'convolution_stepwidth': (0.2, float),
		'convolution_kernelwidth': (2, float),
		'convolution_derivative': (0, int),
		'convolution_amplitude': (1, float),
		'convolution_function': ('Off', str),
		'convolution_widthgauss': (1, float),
		'convolution_widthlorentz': (1, float),

		'series_currenttab': (0, int),
		'series_references': ([], list),
		'series_blendwidth': (0, float),
		'series_qns': (4, int),
		'series_blendminrelratio': (0, float),

		'fit_fitmethod': ('Pgopher', str),
		'fit_uncertainty': (0.05, float),
		'fit_uncertaintystep': (0.01, float),
		'fit_copytoclipboard': (True, bool),
		'fit_xpoints': (1000, int),
		'fit_peakdirection': (1, int),
		'fit_comment': ('', str),
		'fit_polynomrank': (2, int),
		'fit_polynommaxrank': (10, int),
		'fit_offset': (True, bool),

		'isvisible_matplotlibtoolbar': (False, bool),
		'isvisible_controlsmainplot': (True, bool),
		'isvisible_controlsrowscols': (True, bool),
		'isvisible_controlswidth': (True, bool),
		
		'color_exp': ('#ffffff', Color),
		'color_cat': ('#785ef0', Color),
		'color_lin': ('#648fff', Color),
		'color_ref': ('#dc267f', Color),
		'color_fit': ('#fe6100', Color),

		'flag_expformats': ({}, dict),
		'flag_catformats': ({}, dict),
		'flag_linformats': ({}, dict),

		'flag_saveformat':  ({}, dict),

		'flag_tableformatint': ('.0f', str),
		'flag_tableformatfloat': ('.2f', str),

		'flag_extensions': ({"exp": [".csv"], "cat": [".cat"], "lin": [".lin"], "project": [".files"]}, dict),
		'flag_notificationtime': (2000, int),
		'flag_logmaxrows': (1000, int),
		'flag_statusbarmaxcharacters': (100, int),
		'flag_xformatfloat': (".4f", str),
		'flag_syncreferencestocolumns': (True, bool),
		'flag_appendonsave': (True, bool),
		'flag_showmainplotposition': (True, bool),
		'flag_pyckettquanta': (6, int),
		'flag_showseriesarrows': (True, bool),
		'flag_keeponlylastassignment': (False, bool),
		'flag_autoreloadfiles': (True, bool),

		'commandlinedialog_commands': ([], list),
		'commandlinedialog_current': (0, int),
	
		'closebylines_catfstring': ('{x:12.4f} {qns} {ylog}', str),
		'closebylines_linfstring': ('{x:12.4f} {qns}', str),

		'residuals_defaultcolor': ("#000000", Color),
		'residuals_query': ("", str),
		'residuals_colorinput': ("", str),
		'residuals_xvariable': ("", str),
		'residuals_yvariable': ("", str),
		'residuals_autoscale': (True, bool),
		'residuals_blends': (False, bool),

		'blendedlines_lineshape': ("Gauss", str),
		'blendedlines_derivative': (0, int),
		'blendedlines_transparency': (0.2, float),
		'blendedlines_maxfwhm': (10, float),
		'blendedlines_polynom': (0, int),
		'blendedlines_fixedwidth': (False, bool),
		'blendedlines_showbaseline': (True, bool),
		'blendedlines_xpoints': (1000, int),
		'blendedlines_color_total': ("#3d5dff", Color),
		'blendedlines_color_points': ("#ff3352", Color),
		'blendedlines_color_baseline': ("#f6fa14", Color),
		'blendedlines_autopositionpeaks': (True, bool),

		'report_blends': (True, bool),
		'report_query': ('', str),
		
		'seriesfinder_start': ("", str),
		'seriesfinder_stop': ("", str),
		'seriesfinder_results': (10, int),
		'seriesfinder_condition': ("", str),
		'seriesfinder_atype': (True, bool),
		'seriesfinder_btype': (True, bool),
		'seriesfinder_ctype': (True, bool),
		'seriesfinder_xtype': (True, bool),
		'seriesfinder_onlyunassigned': (True, bool),
		
		'peakfinder_peakcolor': ("#4287f5", Color),
		'peakfinder_kwargs': ({}, dict),
		'peakfinder_onlyunassigned': (True, bool),
		'peakfinder_width': (1, float),
		'peakfinder_maxentries': (1000, int),
		
		'energylevels_defaultcolor': ("#000000", Color),
		'energylevels_query': ("", str),
		'energylevels_colorinput': ("", str),
		'energylevels_xvariable': ("", str),
		'energylevels_yvariable': ("", str),
		'energylevels_autoscale': (True, bool),
		
		'cmd_current': (0, int),
		'cmd_commands': ([], list),

		'assignall_fitwidth': (4, float),
		'assignall_maxfwhm': (1, float),
		'assignall_maxdegree': (3, int),

		
		'asap_query': ('', str),
		'asap_resolution': (6e-6, float),
		'asap_weighted': (True, bool),
		'asap_catunitconversionfactor': (0, float),
		'asap_detailviewerwidth': (0, float),
		'asap_detailviewerfilter': (False, bool),
		'asap_assigntransitions': (True, bool),
	}

	def __init__(self, signal, *args, **kwargs):
		super().__init__({key: value[0] for key, value in self.initial_values.items()}, *args, **kwargs)
		self._group_callbacks_counter = AtomicCounter()
		self._grouped_callbacks = set()
		self._grouped_callbacks_lock = threading.Lock()
		self.valuechanged = signal
		self.valuechanged.connect(self.callback)
		self.callbacks = pd.DataFrame(columns=["id", "key", "widget", "function"], dtype="object").astype({"id": np.uint})

	def __setitem__(self, key, value, widget=None):
		if self.get(key) != value:
			super().__setitem__(key, value)
			self.valuechanged.emit((key, value, widget))

	def callback(self, args):
		key, value, widget = args
		if widget:
			callbacks = self.callbacks.query(f"key == @key and widget != @widget")["function"].values
		else:
			callbacks = self.callbacks.query(f"key == @key")["function"].values
		
		counter_value = self._group_callbacks_counter.get_value()
		if counter_value:
			with self._grouped_callbacks_lock:
				self._grouped_callbacks.update(callbacks)
		else:
			for callback in callbacks:
				callback()

	def register(self, keys, function):
		if not isinstance(keys, (tuple, list)):
			keys = [keys]
		for key in keys:
			# id is only needed for callback tied to widgets
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

	def load(self):
		fname = llwpfile(".ini")
		config_parser = configparser.ConfigParser(interpolation=None)
		config_parser.read(fname)

		self.messages = []
		for section in config_parser.sections():
			for key, value in config_parser.items(section):
				fullkey = f"{section.lower()}_{key.lower()}"
				if fullkey in self.initial_values:
					try:
						class_ = self.initial_values[fullkey][1]
						if class_ in (dict, list, tuple):
							value = json.loads(value)
						elif class_ == bool:
							value = True if value in ["True", "1"] else False
						elif class_ == str:
							value = value.encode("utf-8").decode("unicode_escape")
						value = class_(value)
						self[fullkey] = value
					except Exception as E:
						message = f"The value for the option {fullkey} from the option file was not understood."
						self.messages.append(message)
						print(message)
				else:
					self[fullkey] = value
		
		# Special case changing colors for better contrast
		for key, value in self.items():
			if key in self.initial_values and self.initial_values[key][1] == Color:
				if is_dark_theme():
					if matplotlib.colors.to_hex(value) == "#000000":
						self[key] = "#ffffff"
						self.messages.append(f"Changed the color of '{key}' from black to white as it is otherwise invisible.")
				else:
					if matplotlib.colors.to_hex(value) == "#ffffff":
						self[key] = "#000000"
						self.messages.append(f"Changed the color of '{key}' from white to black as it is otherwise invisible.")
	
	def save(self):
		output_dict = {}
		for key, value in self.items():
			category, name = key.split("_", 1)
			category = category.capitalize()
			if category not in output_dict:
				output_dict[category] = {}
			if type(value) in (dict, list, tuple):
				value = json.dumps(value)
			elif type(value) == str:
				value = value.encode("unicode_escape").decode("utf-8")

			output_dict[category][name] = value

		config_parser = configparser.ConfigParser(interpolation=None)
		for section in output_dict:
			config_parser.add_section(section)
			for key in output_dict[section]:
				config_parser.set(section, key, str(output_dict[section][key]))

		with open(llwpfile(".ini"), "w+", encoding="utf-8") as file:
			config_parser.write(file)
		notify_info.emit("Options were saved successfully!")

	def __enter__(self):
		self._group_callbacks_counter.increase()

	def __exit__(self, type, value, traceback):
		counter_value = self._group_callbacks_counter.decrease()
		if counter_value != 0:
			return
		
		with self._grouped_callbacks_lock:
			for callback in self._grouped_callbacks:
				callback()
			self._grouped_callbacks = set()
			
class Geometry():
	_data = {}
	
	@classmethod
	def set(cls, key, value):
		cls._data[key] = value
	
	@classmethod
	def get(cls, key, default=None):
		return(cls._data.get(key, default))

	@classmethod
	def load(cls):
		filename = llwpfile('.geometry')
		if not os.path.isfile(filename):
			return
		
		with open(filename, 'rb') as file:
			cls._data = pickle.load(file)
	
	@classmethod
	def save(cls):
		filename = llwpfile('.geometry')
		with open(filename, 'wb') as file:
			pickle.dump(cls._data, file)

	@classmethod
	def save_widget_geometry(cls, widget):
		key = widget.__class__.__name__
		geometry = widget.geometry()
		cls.set(key, geometry)
	
	@classmethod
	def load_widget_geometry(cls, widget):
		key = widget.__class__.__name__
		geometry = cls.get(key)
		if geometry:
			widget.setGeometry(geometry)

class PlotWidget(QWidget):
	request_redraw = pyqtSignal()

	def __init__(self, parent=None):
		super().__init__(parent)

		self.parent = parent
		self.gui()
		self.from_current_plot()
		self.span_selector = matplotlib.widgets.SpanSelector(self.ax, self.on_range, 'horizontal', useblit=True, button=3)
		self.shortcuts()

	def gui(self):
		layout = QVBoxLayout()
		self.setLayout(layout)

		tmp_layout = QHBoxLayout()
		layout.addLayout(tmp_layout)

		for dir in ("in", "out", "left", "right"):
			tmp_layout.addWidget(QQ(QPushButton, text=dir, change=lambda x, dir=dir: self.move_plot(dir)))

		tmp_layout.addStretch()
		tmp_layout.addWidget(QQ(QPushButton, text="From Current Plot", change=lambda x: self.from_current_plot()))

		self.fig = matplotlib.figure.Figure(dpi=config["plot_dpi"])
		self.plot_canvas = FigureCanvas(self.fig)
		layout.addWidget(self.plot_canvas, 6)
		layout.addStretch()

		self.ax = self.fig.subplots()

		self.exp_coll = matplotlib.collections.LineCollection(np.zeros(shape=(0,2,2)), colors=config["color_exp"], capstyle='round')
		self.cat_coll = matplotlib.collections.LineCollection(np.zeros(shape=(0,2,2)), colors=config["color_cat"], capstyle='round')
		self.lin_coll = self.ax.scatter([], [], color=config['color_ref'], marker="*", zorder=100)
		self.ax.add_collection(self.exp_coll)
		self.ax.add_collection(self.cat_coll)

		self.request_redraw.connect(self.plot_canvas.draw_idle)

	def from_current_plot(self):
		tmp_ax = mainwindow.lwpwidget.get_current_ax()
		self.index = (tmp_ax.row_i, tmp_ax.col_i)
		self.xrange = tmp_ax.xrange

		self.update_plot()

	def update_plot(self):
		scaling = config['plot_yscale']
	
		# Exp Data
		exp_df = ExpFile.get_data(xrange=self.xrange, binning=True)
		self.exp_xs = xs = exp_df["x"].to_numpy()
		self.exp_ys = ys = exp_df["y"].to_numpy()
		exp_yrange = (ys.min(), ys.max()) if len(ys) else (-1, 1)
		segs, colors = [], []

		filenames = exp_df['filename']
		unique_filenames = filenames.unique()
		for unique_filename in unique_filenames:
			mask = (filenames == unique_filename)
			tmp_xs, tmp_ys = xs[mask], ys[mask]

			colors.append( exp_df.loc[mask, 'color'].values )
			if ExpFile.ids[unique_filename].is_stickspectrum:
				segs.append( np.array(((tmp_xs, tmp_xs), (np.zeros(tmp_ys.shape), tmp_ys))).T )
			else:
				segs.append( np.array(((tmp_xs[:-1], tmp_xs[1:]), (tmp_ys[:-1], tmp_ys[1:]))).T )

		if segs:
			segs = np.concatenate(segs)
			colors = np.concatenate(colors)

		self.exp_coll.set(segments=segs, colors=colors)

		# Cat Data
		cat_df = CatFile.get_data(xrange=self.xrange, binning=True)
		self.cat_xs = xs = cat_df["x"].to_numpy()
		self.cat_ys = ys = cat_df["y"].to_numpy()
		cat_yrange = (ys.min(), ys.max()) if len(ys) else (-1, 1)

		if scaling == 'Per Plot':
			ys *= exp_yrange[1] / cat_yrange[1]
		elif scaling in ['Global', 'Custom']:
			ys = ys * config['plot_expcat_factor'] * 10 ** config['plot_expcat_exponent']

		if config['convolution_kernel']:
			kernel_ys = np.array(config['convolution_kernel'])
			convolution_padding = len(kernel_ys) // 2
			stepwidth = config['convolution_stepwidth']

			xmin, xmax = self.xrange
			n_xs = (xmax - xmin) // stepwidth + 1
			con_xs = np.arange(-n_xs, n_xs+1) * stepwidth + (xmax + xmin) / 2

			segs, colors = [], []

			filenames = cat_df['filename']
			unique_filenames = filenames.unique()
			for unique_filename in unique_filenames:
				mask = (filenames == unique_filename)
				tmp_xs, tmp_ys = xs[mask], ys[mask]

				hist_ys, _ = np.histogram(tmp_xs, bins=con_xs, weights=tmp_ys)
				con_ys = np.convolve(hist_ys, kernel_ys, 'full')

				con_ys = con_ys[convolution_padding:-convolution_padding]
				con_xs = con_xs[:-1]
				
				colors.append( cat_df.loc[mask, 'color'].values )
				segs.append( np.array(((con_xs[:-1], con_xs[1:]), (con_ys[:-1], con_ys[1:]))).T )

			if segs:
				segs = np.concatenate(segs)
				colors = np.concatenate(colors)

		else:
			segs = np.array(((xs, xs), (ys*0, ys))).T
			colors = cat_df['color'].to_numpy()

		self.cat_coll.set(segments=segs, colors=colors)

		# Lin Data
		lin_df = LinFile.get_data(xrange=self.xrange, binning=True)
		self.lin_xs = xs = lin_df["x"].to_numpy()
		tuples = [(x, 0) for x in xs]
		tuples = tuples if len(tuples)!=0 else [[None,None]]
		colors = lin_df['color'].to_numpy()

		self.lin_coll.set_offsets(tuples)
		self.lin_coll.set_color(colors)

		# Set Limits
		self.ax.set_xlim(self.xrange)

		if scaling == 'Per Plot':
			yrange = exp_yrange
		elif scaling == 'Global':
			yrange = ExpFile.yrange
		else:
			yrange = (config['plot_yscale_min'], config['plot_yscale_max'])
		margin = config['plot_ymargin']
		yrange = [yrange[0]-margin*(yrange[1]-yrange[0]), yrange[1]+margin*(yrange[1]-yrange[0])]
		if np.isnan(yrange[0]) or np.isnan(yrange[1]) or yrange[0] == yrange[1]:
			yrange = [-1,+1]
		self.ax.set_ylim(yrange)

		self.request_redraw.emit()

	def move_plot(self, dir, factor=None):
		xmin, xmax = self.xrange
		width = xmax - xmin
		center = (xmax + xmin)/2

		if dir == "in":
			width /= 2
		elif dir == "out":
			width *= 2
		elif dir == "left":
			center -= width/2
		elif dir == "right":
			center += width/2
		elif dir == "sin":
			width *= 3/4
		elif dir == "sout":
			width /= 3/4
		elif dir == "sleft":
			center -= width/4
		elif dir == "sright":
			center += width/4
		elif dir == "wheel" and factor != None:
			width *= factor

		self.xrange = center - width/2, center + width/2
		self.update_plot()

	def change_index(self, pos, dir):
		index = list(self.index)
		index[pos] += dir
		
		shape = mainwindow.lwpwidget.lwpaxes.shape
		is_valid_plot = (0 <= index[0] < shape[0] and 0 <= index[1] < shape[1])

		if is_valid_plot:
			self.index = tuple(index)
			tmp_ax = mainwindow.lwpwidget.lwpaxes[self.index]
			xmin, xmax = self.xrange
			curr_width = xmax - xmin

			xmin, xmax = tmp_ax.xrange
			new_center = (xmax + xmin)/2

			self.xrange = (new_center - curr_width/2, new_center + curr_width/2)
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
			"Alt+w": lambda: self.change_index(0, +1),
			"Alt+s": lambda: self.change_index(0, -1),
			"Alt+a": lambda: self.change_index(1, -1),
			"Alt+d": lambda: self.change_index(1, +1),
		}

		for key, function in shortcuts_dict.items():
			tmp = QShortcut(key, self.parent)
			tmp.setContext(Qt.ShortcutContext.WidgetWithChildrenShortcut)
			tmp.activated.connect(function)

	def on_range(self, xmin, xmax):
		xmin_ax, xmax_ax = self.xrange
		if xmax == xmin or xmax > xmax_ax or xmin < xmin_ax:
			return
		
		self.xrange = (xmin, xmax)
		self.update_plot()


drawplot_decorator = DynamicDecorator()
matplotlib_lock = threading.RLock()

##
## Customized PyQt classes
##
class QThread(QThread):
	active_runnables = {}
	lock = threading.RLock()
	threads = set()

	@classmethod
	def threaded_d(cls, func):
		def wrapper(*args, **kwargs):
			thread = cls(func, *args, **kwargs)
			thread.start()
			return(thread)
		return(wrapper)

	def __init__(self, function, *args, **kwargs):
		super().__init__()
		self.function, self.args, self.kwargs = function, args, kwargs
		self.threads.add(self)
		self.finished.connect(lambda x=self: self.threads.remove(x))

	def earlyreturn(self):
		if self.active_runnables[self.function] != self.thread_id:
			raise EarlyReturnError

	def run(self):
		try:
			with self.lock:
				self.thread_id = threading.current_thread().ident
				self.active_runnables[self.function] = self.thread_id
			self.function(*self.args, **self.kwargs, thread=self)
		except EarlyReturnError as E:
			pass
		except Exception as E:
			raise

class QTableWidget(QTableWidget):
	def keyPressEvent(self, event):
		if not csv_copypaste(self, event):
			super().keyPressEvent(event)

class QTableView(QTableView):
	def keyPressEvent(self, event):
		if not csv_copypaste(self, event):
			super().keyPressEvent(event)

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

class QDoubleSpinBoxFullPrec(QDoubleSpinBox):
	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		self.setDecimals(100)
		# AdaptiveDecimalStepType is not implemented in earlier versions of PyQt5
		try:
			self.setStepType(QAbstractSpinBox.StepType.AdaptiveDecimalStepType)
		except:
			pass

	def textFromValue(self, value):
		return(str(value))

class QHBoxLayout(QHBoxLayout):
	def __init__(self, *args, **kwargs):
		margin = kwargs.pop("margin", False)
		super().__init__(*args, **kwargs)
		if not margin:
			self.setContentsMargins(0, 0, 0, 0)

class QVBoxLayout(QVBoxLayout):
	def __init__(self, *args, **kwargs):
		margin = kwargs.pop("margin", False)
		super().__init__(*args, **kwargs)
		if not margin:
			self.setContentsMargins(0, 0, 0, 0)

class QGridLayout(QGridLayout):
	def __init__(self, *args, **kwargs):
		margin = kwargs.pop("margin", False)
		super().__init__(*args, **kwargs)
		if not margin:
			self.setContentsMargins(0, 0, 0, 0)

class FigureCanvas(FigureCanvas):
	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		self.wheelEvent = lambda event: event.ignore()
		# self.setStyleSheet('background-color: #00000000')

		app = QApplication.instance()
		app.styleHints().colorSchemeChanged.connect(self.update_theme)
		self.update_theme()

	def update_theme(self):
		background = 'black' if is_dark_theme() else 'white'
		textcolor = 'white' if is_dark_theme() else 'black'
		
		figure = self.figure

		for ax in figure.get_axes():
			ax.tick_params(color=textcolor, labelcolor=textcolor)
			for spine in ax.spines.values():
				spine.set_edgecolor(textcolor)

		self.setStyleSheet(f'background-color: {background}')
		self.draw_idle()


# class EQDockWidgetMeta(type(QDockWidget)):
	# def __call__(cls, *args, **kwargs):
		# if hasattr(cls, '_instance'):
			# cls._instance.show()
			# cls._instance.raise_()
			# cls._instance.activateWindow()
			# return(cls._instance)
		
		# obj = cls.__new__(cls, *args, **kwargs)
		# cls._instance = obj
		# cls.__init__(obj, *args, **kwargs)
		
		# return(obj)

# class EQDockWidget(QDockWidget, metaclass=EQDockWidgetMeta):
	# default_position = 2
	# default_visible = False


class EQDockWidget(QDockWidget):
	default_position = 2
	default_visible = False

	def __init__(self, *args, **kwargs):
		super().__init__(mainwindow, *args, **kwargs)
		Geometry.load_widget_geometry(self)
		self.is_instantiated = True
		self.show()
		self.activateWindow()
		self.raise_()
		self.setObjectName(self.__class__.__name__)

		if self.default_position is None:
			self.setFloating(True)
		else:
			mainwindow.addDockWidget(Qt.DockWidgetArea(self.default_position), self)
		self.setVisible(self.default_visible)
		
		tmp = QShortcut("Esc", self)
		tmp.setContext(Qt.ShortcutContext.WidgetWithChildrenShortcut)
		tmp.activated.connect(self.close)
		self.__class__.instance = self

	def moveEvent(self, *args, **kwargs):
		Geometry.save_widget_geometry(self)
		return super().moveEvent(*args, **kwargs)

	def resizeEvent(self, *args, **kwargs):
		Geometry.save_widget_geometry(self)
		return super().resizeEvent(*args, **kwargs)

	def show(self, *args, **kwargs):
		screen_box = self.screen().geometry()
		widget_top_left = self.geometry().topLeft()
		widget_bottom_right = self.geometry().bottomRight()
		
		if not (screen_box.contains(widget_top_left) and screen_box.contains(widget_bottom_right)):
			primary_screen = QApplication.instance().primaryScreen()
			self.move(primary_screen.geometry().center()- self.rect().center())
		
		return(super().show(*args, **kwargs))


class QDialog(QDialog):
	def __init__(self, *args, **kwargs):
		super().__init__(mainwindow)
		Geometry.load_widget_geometry(self)
		QShortcut("Esc", self).activated.connect(lambda: self.done(0))

	def moveEvent(self, *args, **kwargs):
		Geometry.save_widget_geometry(self)
		return super().moveEvent(*args, **kwargs)

	def resizeEvent(self, *args, **kwargs):
		Geometry.save_widget_geometry(self)
		return super().resizeEvent(*args, **kwargs)


##
## File claasses
##

class SpecialFilesHandler():
	@classmethod
	def save_files(cls, dict_):
		return(dict_)
	
	@classmethod
	def load_files(cls, dict_):
		return(dict_)
	
	@classmethod
	def sort_files_by_type(cls, files):
		return(files)

class File():
	ids = {}
	lock = threading.RLock()

	default_color_key = ""
	has_y_data = True

	yrange = [0, 0]
	special_file_handler = SpecialFilesHandler()

	additional_dtypes = {
		'x0': np.float64,
		'y0': np.float64,
		'color': str,
		'visible': bool,
		'filename': 'category',
	}
	files_watcher = QFileSystemWatcher()
	files_watcher.fileChanged.connect(lambda x: File.autoreload_file(x))

	@classmethod
	def autoreload_file(cls, filename):
		if not config['flag_autoreloadfiles']:
			return

		sorted_files = cls.sort_files_by_type([filename])
		for filetype, files in sorted_files.items():
			if not files:
				continue

			class_ = {'exp': ExpFile, 'cat': CatFile, 'lin': LinFile}.get(filetype)
			if not class_:
				return

			file = class_.ids.get(filename)
			if not file:
				return
			
			if os.path.getsize(filename) == 0:
				return
			
			file.load_file()

	def __new__(cls, filename, *args, load_manually=False, **kwargs):
		if cls == NewAssignments:
			id = filename
		else:
			id = os.path.abspath(filename)
		if id in cls.ids.keys():
			existing_file = cls.ids[id]
			if not load_manually:
				existing_file.load_file()
			return(existing_file)
		
		return super().__new__(cls)

	def __init__(self, filename, load_manually=False, default_values={}):
		if hasattr(self, 'is_initialized'):
			return
		
		filename = os.path.abspath(filename)
		self.filename_abs = filename
		self.dirname_abs, self.basename = os.path.split(filename)
		self.extension = os.path.splitext(filename)[1]
		
		if self.filename_abs not in self.files_watcher.files():
			self.files_watcher.addPath(self.filename_abs)

		self.ids[self.filename_abs] = self

		self.gui_widgets = {}
		self.more_settings_dialog = None
		self._color = None
		self._is_visible = None
		self.set_default_values(default_values)
		
		if not load_manually:
			self.load_file()
		
		if hasattr(FileWindow, 'instance'):
			FileWindow.instance.fileaddition_requested.emit(self.__class__, self.filename_abs)
		self.is_initialized = True

	def set_default_values(self, default_values={}):
		self.color = default_values.get('color', config[self.default_color_key])
		self.color_query = default_values.get('color_query')
		self.query = default_values.get('query')
		self.is_visible = default_values.get('visible', True)
		self.xtransformation = default_values.get('xtransformation')
		self.ytransformation = default_values.get('ytransformation')
		
		if hasattr(self, 'is_stickspectrum'):
			self.is_stickspectrum = default_values.get('is_stickspectrum', False)

	def apply_all(self):
		self.apply_color()
		self.apply_visibility()
		self.apply_xtransformation()
		if self.has_y_data:
			self.apply_ytransformation()

	def check_file(self):
		fname = self.filename_abs
		if not os.path.isfile(fname):
			notify_error.emit(f"The file {fname} could not be found. Please check the file.")
			raise CustomError(f"The file {fname} could not be found. Please check the file.")
		
		if os.path.getsize(fname) == 0:
			notify_warning.emit(f"The file {fname} is empty and was therefore skipped.")
			raise CustomError(f"The file {fname} is empty and was therefore skipped.")
	
	@QThread.threaded_d
	@status_d
	def load_file(self, thread=None):
		with self.sort_df_counter:

			self.check_file()
			data = self.load_file_core()

			# Add additional columns
			data['x0'] = data['x']
			if self.has_y_data:
				data['y0'] = data['y']
			data['color'] = self.color
			data['visible'] = True
			
			with self.lock:
				df = self.get_df()
				# Important to explicitly convert filename to category again; Otherwise it will be an object which makes masking 
				# a factor of 1000 slower
				self.set_df(pd.concat([df[df['filename'] != self.filename_abs], data], ignore_index=True).astype({'filename': 'category'}))
				self.clean_up_data()
				self.apply_all()

			# df = self.get_df()
			# mask = (df['filename'] == self.filename_abs)
			# xs = df.loc[mask, 'x']
			# self.xmin, self.xmax = xs.min(), xs.max()

			self.clear_caches()

			if not isinstance(self, NewAssignments):
				notify_info.emit(f"Successfully loaded '{self.filename_abs}'")
	
	@classmethod
	def clean_up_data(cls):
		pass
	
	@classmethod
	@lru_cache
	def has_results(cls, query):
		with cls.lock:
			n_results = len(cls.get_data().query(query))
		return(n_results)
	
	@classmethod
	@lru_cache
	def query_c(cls, query):
		with cls.lock:
			return(cls.get_data().query(query))

	@classmethod
	def clear_caches(cls):
		cls.query_c.cache_clear()
		cls.has_results.cache_clear()

	def to_dict(self):
		dict_ = {
			'filename': self.filename_abs,
			'color': self.color,
			'color_query': self.color_query,
			'query': self.query,
			'visible': self.is_visible,
			'xtransformation': self.xtransformation,
			'ytransformation': self.ytransformation,
		}
		if hasattr(self, 'is_stickspectrum'):
			dict_['is_stickspectrum'] = self.is_stickspectrum
		return(dict_)

	@classmethod
	def save_files_gui(cls):
		savepath, filter = QFileDialog.getSaveFileName(None, 'Project file', '')
		if not savepath:
			return

		cls.save_files(savepath)
		notify_info.emit(f'Saved the current files as project \'{savepath}\'.')

	@classmethod
	def save_files(cls, filename):
		dict_ = {}
		for subclass in cls.__subclasses__():
			tmp = [instance.to_dict() for instance in subclass.ids.values() if not isinstance(instance, NewAssignments)] 
			dict_[subclass.__name__] = tmp
		
		dict_ = cls.special_file_handler.save_files(dict_)

		with open(filename, 'w+') as file:
			json.dump(dict_, file, indent=4)

	@classmethod
	def load_files_gui(cls):
		filename, filter = QFileDialog.getOpenFileName(None, 'Choose Project to load')
		if not filename:
			return
		
		thread = cls.load_files(filename)
		thread.wait()
		notify_info.emit(f'Loaded the project \'{filename}\'.')

	@classmethod
	@QThread.threaded_d
	def load_files(cls, filename, thread=None):
		with open(filename, 'r') as file:
			dict_ = json.load(file)
		
		dict_ = cls.special_file_handler.load_files(dict_)

		label_to_class = {subclass.__name__: subclass for subclass in cls.__subclasses__()}
		threads = []
		for label, filedicts in dict_.items():
			subclass = label_to_class[label]
			for filedict in filedicts:
				filename = filedict.pop('filename')
				tmp_file = subclass(filename, default_values=filedict, load_manually=True)
				threads.append(tmp_file.load_file())

		for thread_ in threads:
			thread_.wait()
		
	@classmethod
	def get_data(cls, xrange=None, binning=False):
		with cls.lock:
			df = cls.df

			if xrange is not None:
				x_start = df["x"].searchsorted(xrange[0], side="left")
				x_stop  = df["x"].searchsorted(xrange[1], side="right")
				df = df.iloc[x_start:x_stop].copy()

			df = df[df['visible']]

			if binning and xrange is not None:
				bins = config['plot_bins']
				nobinning = config['plot_skipbinning']
				binwidth = (xrange[1]-xrange[0]) / bins

				if len(df) > max(bins, nobinning)  and binwidth != 0:
					df = bin_data(df, binwidth, xrange)

		return(df)

	@classmethod
	def get_df(cls):
		return(cls.df)
	
	@classmethod
	def set_df(cls, value):
		cls.df = value

	@classmethod
	def xs_to_indices(cls, xmins, xmaxs):
		with cls.lock:
			minindices = cls.df["x"].searchsorted(xmins, side="left")
			maxindices = cls.df["x"].searchsorted(xmaxs, side="right")
		return(minindices, maxindices)

	@classmethod
	def sort_df(cls):
		with cls.lock:
			cls.df = cls.df.sort_values('x')
			if cls.has_y_data:
				tmp = cls.df[ cls.df['visible'] ]['y']
				cls.yrange = tmp.min(), tmp.max()
		mainwindow.lwpwidget.set_data()

	@classmethod
	def sort_files_by_type(cls, files):
		files = cls.special_file_handler.sort_files_by_type(files)
		types = config['flag_extensions']
		files_by_type = {key: [] for key in list(types.keys())}

		for file in files:
			path, extension = os.path.splitext(file)
			extension = extension if extension else os.path.basename(path)
			type = None
			for key, value in types.items():
				if extension in value:
					type = key
					break

			if type is None:
				item, ok = QInputDialog.getItem(None, "Choose File Type", f"Choose the file type for the extension \"{extension}\":", [x.capitalize() for x in types], editable=False)
				if not (ok and item):
					continue
				types[item.lower()].append(extension)
				type = item.lower()

			files_by_type[type].append(file)
		return(files_by_type)

	@classmethod
	@QThread.threaded_d
	@status_d
	def add_multiple(cls, filepaths, thread=None):
		if not filepaths:
			return
		
		with cls.sort_df_counter:
			files = [cls(file, load_manually=True) for file in filepaths]
			dataframes = []
			filenames = []
			exceptions = {}
			
			for file in files:
				try:
					dataframes.append(file.load_file_core())
					filenames.append(file.filename_abs)
				except Exception as E:
					raise E
					exceptions[file.filename_abs] = str(E)
			
			data = pd.concat(dataframes)
			data['x0'] = data['x']
			if cls.has_y_data:
				data['y0'] = data['y']
			data['color'] = default_color = config[cls.default_color_key]
			data['visible'] = True

			
			with cls.lock:
				df = cls.get_df()
				# Important to explicitly convert filename to category again; Otherwise it will be an object which makes masking 
				# a factor of 1000 slower
				cls.set_df(pd.concat([df[~df['filename'].isin(filenames)], data], ignore_index=True).astype({'filename': 'category'}))
				cls.clean_up_data()
				
				# Apply color
				colors = [file.color for file in files]
				color_queries = [file.color_query for file in files]
				
				for file in files:
					if not (file.color == default_color and not file.color_query):
						file.apply_color()
					
					if not (file.is_visible and not file.query):
						file.apply_visibility()
					
					if file.xtransformation:
						file.apply_xtransformation()
					
					if file.has_y_data and file.ytransformation:
						file.apply_ytransformation()

			cls.clear_caches()
		
		if len(exceptions) == 0:
			notify_info.emit(f"Successfully loaded all files.")
		else:
			error_string = '\n'.join(f'{file}: {error}' for file, error in exceptions.items())
			notify_warning.emit(f"There were errors when loading the requested files:\n{error_string}")


	@staticmethod
	def add_files_dialog():
		files = QFileDialog.getOpenFileNames(None, f'Choose File(s)',)[0]
		files_by_type = File.sort_files_by_type(files)
		File.add_multiple_files_by_type(files_by_type)

	@staticmethod
	@QThread.threaded_d
	def add_multiple_files_by_type(files_by_type, thread=None):
		threads = (
			ExpFile.add_multiple(files_by_type.get("exp", [])),
			CatFile.add_multiple(files_by_type.get("cat", [])),
			LinFile.add_multiple(files_by_type.get("lin", [])),
		)

		for thread_ in threads:
			thread_.wait()

		threads = (File.load_files(project) for project in files_by_type.get('project', []))
		for thread_ in threads:
			thread_.wait()
		
		CatFile.check_series_qns()

	@classmethod
	def gui_settings_general(cls):
		layout = QVBoxLayout()
		

		buttons_dict = {
			'Reread': cls.reread_all,
			'Reset': cls.reset_all,
			'Hide/Show': cls.toggle_all,
			'Delete': cls.delete_all,
		}

		dcolor = config[cls.default_color_key]
		stylesheet = f"background-color: {Color.rgbt_to_trgb(dcolor)}"
		color_label = QQ(QLabel, text='Default Color: ')
		color_input = QQ(QLineEdit, text=dcolor, maxWidth=200, change=cls.gui_change_color_all)
		color_picker = QQ(QToolButton, stylesheet=stylesheet, change=cls.gui_change_color_all)
		
		tmp_layout = QHBoxLayout()
		layout.addLayout(tmp_layout)
		
		for key, callback in buttons_dict.items():
			tmp_layout.addWidget(QQ(QPushButton, text=key, change=lambda x, callback=callback: callback()))
		tmp_layout.addStretch(1)
		
		tmp_layout = QHBoxLayout()
		layout.addLayout(tmp_layout)
		
		for widget in (color_label, color_input, color_picker):
			tmp_layout.addWidget(widget)
		tmp_layout.addStretch(1)
		
		cls.gui_widgets_class = {
			'color_picker': color_picker,
			'color_input': color_input,
		}

		return(layout)

	def gui_settings_widgets(self):
		dcolor = self.color
		stylesheet = f"background-color: {Color.rgbt_to_trgb(dcolor)}"
		toggle_text = 'Hide' if self.is_visible else 'Show'
		
		widgets = {
			'label': QQ(QLabel, text=self.basename, tooltip=self.filename_abs, enabled=self.is_visible),
			'color_input': QQ(QLineEdit, text=dcolor, maxWidth=200, change=self.gui_change_color),
			'color_picker': QQ(QToolButton, stylesheet=stylesheet, change=self.gui_change_color),
			'toggle_visibility': QQ(QToolButton, text=toggle_text, change=self.gui_toggle_visbility),
			'settings_dialog': QQ(QToolButton, text='⚙', change=self.gui_more_settings_dialog),
			'reread': QQ(QToolButton, text='⟲', change=lambda x: self.load_file()),
			'delete': QQ(QToolButton, text='×', change=lambda x: self.gui_delete()),
		}

		for widget in self.gui_widgets.values():
			widget.setParent(None)

		self.gui_widgets = widgets
		return(widgets.values())
		
	def gui_more_settings_dialog(self):
		self.more_settings_dialog = FileAdditionalSettingsDialog.show_dialoag(self)
			
	@classmethod
	def gui_change_color_all(cls, argument):
		dcolor = config[cls.default_color_key]
		if type(argument) == str:
			color = argument
		else:
			color = QColorDialog.getColor(initial=QColor(Color.rgbt_to_trgb(dcolor)), options=QColorDialog.ColorDialogOption.ShowAlphaChannel)
			if color.isValid():
				color = Color.trgb_to_rgbt(color.name(QColor.NameFormat.HexArgb))
			else:
				return
		
		try:
			color = Color(color)
		except CustomError:
			return

		config[cls.default_color_key] = color

		color_input = cls.gui_widgets_class.get('color_input')
		if color_input and color_input.text() != color:
			color_input.setText(color)
		color_picker = cls.gui_widgets_class.get('color_picker')
		if color_picker:
			color_picker.setStyleSheet(f"background-color: {Color.rgbt_to_trgb(color)}")
		
	def gui_change_color(self, argument):
		if type(argument) == str:
			color = argument
		else:
			color = QColorDialog.getColor(initial=QColor(Color.rgbt_to_trgb(self.color)), options=QColorDialog.ColorDialogOption.ShowAlphaChannel)
			if color.isValid():
				color = Color.trgb_to_rgbt(color.name(QColor.NameFormat.HexArgb))
			else:
				return
		
		try:
			color = Color(color)
		except CustomError:
			return
		
		self.color = color
		self.apply_color()
		mainwindow.lwpwidget.set_data()

	def gui_toggle_visbility(self, _):
		self.toggle_visibility()
		self.apply_visibility()
		mainwindow.lwpwidget.set_data()

	@classmethod
	def reread_all(cls):
		if cls == File:
			for subcls in cls.__subclasses__():
				subcls.reread_all()
		else:
			for file in cls.ids.values():
				file.load_file()
	
	@classmethod
	def reset_all(cls):
		for file in cls.ids.values():
			file.set_default_values()
		
		df = cls.df
		df['color'] = config[cls.default_color_key]
		df['visible'] = True
		df['x'] = df['x0']
		if cls.has_y_data:
			df['y'] = df['y0']
				
		mainwindow.lwpwidget.set_data()

	
	@classmethod
	def delete_all(cls):
		for file in list(cls.ids.values()):
			file.delete()
		mainwindow.lwpwidget.set_data()
	
	@classmethod
	def toggle_all(cls):
		files = cls.ids.values()
		all_are_visible = all([file.is_visible for file in files])
		visibility = not all_are_visible

		for file in cls.ids.values():
			file.is_visible = visibility

		if visibility:
			for file in cls.ids.values():
				file.apply_visibility()
		else:
			df = cls.get_df()
			df['visible'] = False
			cls.set_df(df)
		
		mainwindow.lwpwidget.set_data()


	@property
	def color(self):
		return(self._color)
	
	@color.setter
	def color(self, value):
		self._color = value
		
		color_input = self.gui_widgets.get('color_input')
		if color_input and color_input.text() != value:
			color_input.setText(value)
		color_picker = self.gui_widgets.get('color_picker')
		if color_picker:
			color_picker.setStyleSheet(f"background-color: {Color.rgbt_to_trgb(value)}")

	def apply_color(self):
		df = self.__class__.df
		mask = (df['filename'] == self.filename_abs)

		df.loc[mask, 'color'] = self.color
		
		# Color query
		if not self.color_query:
			self.clear_caches()
			return
		
		for command in self.color_query.split("\n"):
			if not command.strip():
				continue
			color, query = command.split(";")
			indices = (df.loc[mask].query(query)).index
			df.loc[indices, "color"] = color
		self.clear_caches()


	def toggle_visibility(self):
		self.is_visible = not self.is_visible

	@property
	def is_visible(self):
		return(self._is_visible)
	
	@is_visible.setter
	def is_visible(self, value):
		self._is_visible = value

		widget = self.gui_widgets.get('toggle_visibility')
		if widget:
			widget.setText('Hide' if self.is_visible else 'Show')
		label = self.gui_widgets.get('label')
		if label:
			label.setEnabled(self.is_visible)

	def apply_visibility(self):
		df = self.__class__.df
		mask = (df['filename'] == self.filename_abs)

		if self.is_visible and not self.query:
			df.loc[mask, 'visible'] = True
		elif self.is_visible and self.query:
			df.loc[mask, 'visible'] = df.loc[mask].eval(self.query)
		else:
			df.loc[mask, 'visible'] = False
		self.clear_caches()
		
	def apply_transformation(self, col, transform, fallback_col):
		df = self.__class__.df
		mask = (df['filename'] == self.filename_abs)
		
		if not transform:
			df.loc[mask, col] = df.loc[mask, fallback_col]
		else:
			df.loc[mask, col] = df.loc[mask].eval(transform)
		self.clear_caches()


	def apply_xtransformation(self, *args, **kwargs):
		self.apply_transformation('x', self.xtransformation, 'x0', *args, **kwargs)

	def apply_ytransformation(self, *args, **kwargs):
		self.apply_transformation('y', self.ytransformation, 'y0', *args, **kwargs)

	def gui_delete(self):
		self.delete()
		mainwindow.lwpwidget.set_data()

	def delete(self):
		# Delete from Filewatcher
		resp = self.files_watcher.removePath(self.filename_abs)

		# Delete row from files window
		for widget in self.gui_widgets.values():
			widget.setParent(None)
		self.gui_widgets = None
		
		# Close/delete file
		dialog = self.more_settings_dialog
		if dialog and dialog.isVisible():
			dialog.done(0)

		# Delete data from class dataframe
		with self.lock:
			df = self.__class__.df
			mask = df[df['filename'] == self.filename_abs].index
			self.__class__.df = df.drop(mask)
			del self.ids[self.filename_abs]
		
		self.clear_caches()
		del self
	
class CatFile(File):
	ids = {}
	lock = threading.RLock()
	sort_df_counter = AtomicCounter()

	default_color_key = "color_cat"
	dtypes = {**pyckett.cat_dtypes_from_quanta(N_QNS), **File.additional_dtypes}
	df = pd.DataFrame(columns=dtypes.keys()).astype(dtypes)

	def load_file_core(self):
		if self.extension not in config['flag_catformats']:
			data = pyckett.cat_to_df(self.filename_abs, sort=False)
		else:
			kwargs = config['flag_catformats'][self.extension].copy()
			y_is_log = format.pop('y_is_log', False)

			data = pd.read_fwf(self.filename_abs, **kwargs)
			data['filename'] = self.filename_abs

			if y_is_log:
				data['y'] = 10 ** data['y']
				
			for column in self.dtypes.keys():
				if column.startswith('qn'):
					continue
				
				if column not in data.columns:
					message = f'The format used to load the file \'{self.filename_abs}\' does not provide the column \'{column}\'.'
					notify_warning.emit(message)
					raise GUIAbortedError(message)
			data = data[self.dtypes.keys()]

		return(data)

	@classmethod
	def clean_up_data(cls):
		qn_columns = pyckett.qnlabels_from_quanta(N_QNS)
		cls.df.loc[:, qn_columns] = cls.df.loc[:, qn_columns].fillna(pyckett.SENTINEL).astype(pyckett.pickett_int)
		# self.__class__.df = self.__class__.df.dropna()

	@classmethod
	def check_series_qns(cls):
		df = cls.df
		
		if not len(df):
			return
		
		qnu_labels = [f'qnu{i+1}' for i in range(N_QNS)]
		noq = len(qnu_labels)
		for i, qnu_label in enumerate(qnu_labels):
			unique_values = df[qnu_label].unique()
			if len(unique_values) == 1 and unique_values[0] == pyckett.SENTINEL:
				noq = i
				break
		
		config['series_qns'] = noq
		notify_info.emit(f'After analysing your cat files the number of QNs was set to {noq}.')

CatFile.sort_df_counter.callback = CatFile.sort_df

class LinFile(File):
	ids = {}
	lock = threading.RLock()
	sort_df_counter = AtomicCounter()

	default_color_key = "color_lin"
	dtypes = {**pyckett.lin_dtypes_from_quanta(N_QNS), **File.additional_dtypes}
	del dtypes['y0']
	df = pd.DataFrame(columns=dtypes.keys()).astype(dtypes)
	
	has_y_data = False

	def load_file_core(self):
		if self.extension not in config['flag_linformats']:
			data = pyckett.lin_to_df(self.filename_abs, sort=False)
		else:
			kwargs = config['flag_linformats'][self.extension].copy()

			data = pd.read_fwf(self.filename_abs, **kwargs)
			data['filename'] = self.filename_abs

			for column in self.dtypes.keys():
				if column.startswith('qn'):
					continue
				
				if column not in data.columns:
					message = f'The format used to load the file \'{self.filename_abs}\' does not provide the column \'{column}\'.'
					notify_warning.emit(message)
					raise GUIAbortedError(message)
			
			data = data[self.dtypes.keys()]

		return(data)
	
	@classmethod
	def clean_up_data(cls):
		qn_columns = pyckett.qnlabels_from_quanta(N_QNS)
		cls.df.loc[:, qn_columns] = cls.df.loc[:, qn_columns].fillna(pyckett.SENTINEL).astype(pyckett.pickett_int)
		# self.__class__.df = self.__class__.df.dropna()

LinFile.sort_df_counter.callback = LinFile.sort_df


class ExpFile(File):
	ids = {}
	lock = threading.RLock()
	sort_df_counter = AtomicCounter()

	default_color_key = "color_exp"
	dtypes = {'x': np.float64, 'y': np.float64, **File.additional_dtypes}
	df = pd.DataFrame(columns=dtypes.keys()).astype(dtypes)
	
	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		self.is_stickspectrum = False

	def load_file_core(self):
		kwargs = config['flag_expformats'].get(self.extension, {})
		data = exp_to_df(self.filename_abs, **kwargs)
		return(data)

ExpFile.sort_df_counter.callback = ExpFile.sort_df

class NewAssignments(LinFile):
	_instance = None
	def __init__(self, filename, load_manually=False):
		if hasattr(self, 'is_initialized'):
			return
		self.is_initialized = True
		self.__class__._instance = self

		dtypes = pyckett.lin_dtypes_from_quanta(N_QNS)
		self.new_assignments_df = pd.DataFrame(columns=dtypes.keys()).astype(dtypes)
		self.filename_abs = filename
		self.basename = "New Assignments"
		self.dirname_abs = self.extension = None
		
		self.ids[self.filename_abs] = self
		self.gui_widgets = {}
		self._color = None
		self._is_visible = None
		self.set_default_values()
		
		if not load_manually:
			self.load_file()
		
		if hasattr(FileWindow, 'instance'):
			FileWindow.instance.fileaddition_requested.emit(self.LinFile, self.filename_abs)

	def check_file(self):
		return
	
	def load_file_core(self, thread=None):
		return(self.new_assignments_df.copy())
	
	@classmethod
	def get_df(cls):
		return(LinFile.df)
	
	@classmethod
	def set_df(cls, value):
		LinFile.df = value

	def gui_settings_widgets(self):
		dcolor = self.color
		stylesheet = f"background-color: {Color.rgbt_to_trgb(dcolor)}"
		toggle_text = 'Hide' if self.is_visible else 'Show'
		
		widgets = {
			'label': QQ(QLabel, text=self.basename, tooltip=self.filename_abs, enabled=self.is_visible),
			'color_input': QQ(QLineEdit, text=dcolor, maxWidth=200, change=self.gui_change_color),
			'color_picker': QQ(QToolButton, stylesheet=stylesheet, change=self.gui_change_color),
			'toggle_visibility': QQ(QToolButton, text=toggle_text, change=self.gui_toggle_visbility),
			'settings_dialog': QQ(QToolButton, text='⚙', change=self.gui_more_settings_dialog),
		}

		for widget in self.gui_widgets.values():
			widget.setParent(None)

		self.gui_widgets = widgets
		return(widgets.values())
	
	def delete(self):
		return
	
	def add_row(self, row_dict):
		rows_dict = {key: [value] for key, value in row_dict.items()}
		self.add_rows(rows_dict)

	def add_rows(self, rows_dict):
		dtypes = pyckett.lin_dtypes_from_quanta(N_QNS)
		new_rows = pd.DataFrame(rows_dict).astype(dtypes)[dtypes.keys()]
		new_rows['filename'] = self.filename_abs

		self.new_assignments_df = pd.concat( (self.new_assignments_df, new_rows) )
		if config['flag_keeponlylastassignment']:
			subset = [f'qn{ul}{i+1}' for ul in 'ul' for i in range(config['series_qns'])]
			self.new_assignments_df = self.new_assignments_df.drop_duplicates(subset=subset, keep='last', ignore_index=True)
		self.new_assignments_df = self.new_assignments_df.reset_index(drop=True)

		self.load_file()

		new_assignments_window = NewAssignmentsWindow.instance
		new_assignments_window.model.update()
		new_assignments_window.model.resize_columns()
		new_assignments_window.scroll_bottom()

	@classmethod
	def get_instance(cls):
		if cls._instance is None:
			return(cls('__newassignments__'))
		else:
			return(cls._instance)
	
	def get_new_assignments_df(self):
		return(self.new_assignments_df)
	
	def save_gui(self):
		append = config['flag_appendonsave']
		format = config['flag_saveformat']
		
	
		if sys.platform == 'darwin':
			options = {"options": QFileDialog.Option.DontConfirmOverwrite | QFileDialog.Option.DontUseNativeDialog} if append else {"options": QFileDialog.Option.DontUseNativeDialog}
		else:
			options = {"options": QFileDialog.Option.DontConfirmOverwrite} if append else {}
		savepath, extension = QFileDialog.getSaveFileName(None, 'Save file', '', **options)
		if not savepath:
			return

		self.save(savepath, append, format)
		notify_info.emit(f"The {len(self.new_assignments_df)} new assignments were saved to the file \'{savepath}\'.")

	def save_backup(self):
		savepath = llwpfile('.lin')
		self.save(savepath, False, None)
	
	def save(self, savepath, append, format):
		handle = "a+" if append else "w+"
		df = self.get_new_assignments_df().copy()

		with open(savepath, handle, encoding="utf-8") as file:
			if format:
				np.savetxt(file, df[format.get("names", df.columns)], delimiter=format.get("delimiter", " "), fmt=format.get("format", '%.18e'))
			else:
				file.write(pyckett.df_to_lin(df))

class FileAdditionalSettingsDialog(QDialog):
	open_dialogs = {}

	@classmethod
	def show_dialoag(cls, file):
		if file in cls.open_dialogs:
			dialog = cls.open_dialogs[file]
			dialog.done(0)
			return
		else:
			dialog = cls(file)
			cls.open_dialogs[file] = dialog
			dialog.show()
			return(dialog)
	
	def __init__(self, file):
		super().__init__()

		self.setModal(False)
		self.setWindowTitle(f'Settings for {file.basename}')
		self.file = file
		self.finished.connect(self.on_exit)
		
		self.layout = QVBoxLayout(margin=True)
		self.setLayout(self.layout)

		ph_colorq = "Enter custom color and query to color specific lines differently. E.g. enter '#ff0000; qn1 < 20' to color all levels with the first quantum number below 20 red."

		self.widgets = {
			'query': {
				'label': QQ(QLabel, text='Filter: '),
				'widget': QQ(QPlainTextEdit, value=file.query),
				'button': QQ(QToolButton, text='Update', change=self.update_query),
			},
			'xtransformation': {
				'label': QQ(QLabel, text='x-Transformation: '),
				'widget': QQ(QPlainTextEdit, value=file.xtransformation),
				'button': QQ(QToolButton, text='Update', change=self.update_xtransformation),
			},
			'ytransformation': {
				'label': QQ(QLabel, text='y-Transformation: '),
				'widget': QQ(QPlainTextEdit, value=file.ytransformation),
				'button': QQ(QToolButton, text='Update', change=self.update_ytransformation),
			} if self.file.has_y_data else {},
			'color_query': {
				'label': QQ(QLabel, text='Color Query: '),
				'widget': QQ(QPlainTextEdit, value=file.color_query, placeholder=ph_colorq),
				'button': QQ(QToolButton, text='Update', change=self.update_color_query),
			},
		}

		for _, tmp in self.widgets.items():
			if not tmp:
				continue
			tmp_layout = QHBoxLayout()
			tmp_layout.addWidget(tmp['label'])
			tmp_layout.addStretch(1)
			tmp_layout.addWidget(tmp['button'])

			self.layout.addLayout(tmp_layout)
			self.layout.addWidget(tmp['widget'])
		
		if hasattr(self.file, 'is_stickspectrum'):
			tmp_widget = QQ(QCheckBox, text='Is Stick Spectrum: ', change=self.update_stickspectrum, value=file.is_stickspectrum)
			self.widgets['is_stickspectrum'] = tmp_widget
			self.layout.addWidget(tmp_widget)

		self.button_box = QDialogButtonBox(QDialogButtonBox.StandardButton.Apply | QDialogButtonBox.StandardButton.Reset | QDialogButtonBox.StandardButton.Close)
		self.button_box.rejected.connect(lambda: self.done(0))

		self.button_box.button(QDialogButtonBox.StandardButton.Reset).clicked.connect(lambda _: self.reset_all())
		self.button_box.button(QDialogButtonBox.StandardButton.Apply).clicked.connect(lambda _: self.apply_all())

		self.layout.addWidget(self.button_box)

	def reset_all(self):
		self.widgets['query']['widget'].setPlainText('')
		self.widgets['xtransformation']['widget'].setPlainText('')
		if self.file.has_y_data:
			self.widgets['ytransformation']['widget'].setPlainText('')
		self.widgets['color_query']['widget'].setPlainText('')

		if hasattr(self.file, 'is_stickspectrum'):
			self.widgets['is_stickspectrum'].setChecked(False)
		
		self.apply_all()

	def apply_all(self):
		self.done(0)
		with self.file.sort_df_counter:
			self.update_query()
			self.update_xtransformation()
			if self.file.has_y_data:
				self.update_ytransformation()
			self.update_color_query()
			if hasattr(self.file, 'is_stickspectrum'):
				self.update_stickspectrum()
		
		mainwindow.lwpwidget.set_data()
		

	def update_query(self, _=None):
		file = self.file
		query = self.widgets['query']['widget'].toPlainText()
		file.query = query
		file.apply_visibility()
		mainwindow.lwpwidget.set_data()
		
	def update_xtransformation(self, _=None):
		with self.file.sort_df_counter:
			file = self.file
			xtransformation = self.widgets['xtransformation']['widget'].toPlainText()
			file.xtransformation = xtransformation
			file.apply_xtransformation()
			
		mainwindow.lwpwidget.set_data()

	def update_ytransformation(self, _=None):
		file = self.file
		ytransformation = self.widgets['ytransformation']['widget'].toPlainText()
		file.ytransformation = ytransformation
		file.apply_ytransformation()
		mainwindow.lwpwidget.set_data()

	def update_color_query(self, _=None):
		file = self.file
		color_query = self.widgets['color_query']['widget'].toPlainText()
		file.color_query = color_query
		file.apply_color()
		mainwindow.lwpwidget.set_data()
	
	def update_stickspectrum(self, _=None):
		file = self.file
		is_stickspectrum = self.widgets['is_stickspectrum'].isChecked()
		file.is_stickspectrum = is_stickspectrum
		mainwindow.lwpwidget.set_data()

	def on_exit(self, _=None):
		del self.__class__.open_dialogs[self.file]


##
## Application and MainWindow
##
class LWPAx():
	fit_vline = None
	fit_curve = None
	fit_methods = ('Pgopher', 'Polynom', 'MultiPolynom', 'Gauss', 'Lorentz',
				'Voigt', 'Gauss 1st Derivative', 'Lorentz 1st Derivative', 'Voigt 1st Derivative',
				'Gauss 2nd Derivative', 'Lorentz 2nd Derivative', 'Voigt 2nd Derivative', )

	def __init__(self, ax, row_i, col_i):
		self.ax = ax
		self.row_i = row_i
		self.col_i = col_i
		
		self.ref_position = None
		self.xrange = (-1, 1)
		self.indices = None
		self.annotation = None
		self.qns = None
		
		with matplotlib_lock:
			self.span =  matplotlib.widgets.SpanSelector(ax, lambda xmin, xmax: self.on_range(xmin, xmax), 'horizontal', useblit=True, button=1)
		
		self.exp_coll = matplotlib.collections.LineCollection(np.zeros(shape=(0,2,2)), colors=config["color_exp"], capstyle='round')
		self.cat_coll = matplotlib.collections.LineCollection(np.zeros(shape=(0,2,2)), colors=config["color_cat"], capstyle='round')
		self.lin_coll = ax.scatter([], [], color=config['color_ref'], marker="*", zorder=100)

		with matplotlib_lock:
			ax.add_collection(self.exp_coll)
			ax.add_collection(self.cat_coll)
			
			ax.yaxis.set_visible(False)
			if row_i:
				ax.xaxis.set_visible(False)
			else:
				ax.set_xticks([])

	# @status_d
	@QThread.threaded_d
	@drawplot_decorator.d
	def update(self, thread=None):
		ax = self.ax
		ax.set_xlim(self.xrange)

		self.set_xticklabels()

		yrange_exp = [-1, 1]

		for datatype, cls in {'exp': ExpFile, 'cat': CatFile, 'lin': LinFile}.items():
			minindex, maxindex = self.indices[datatype]
			xmin, xmax = self.xrange
			with cls.lock:
				dataframe = cls.df.iloc[minindex:maxindex].copy()
			bins = config['plot_bins']
			nobinning = config['plot_skipbinning']
			scaling = config["plot_yscale"]

			dataframe = dataframe[dataframe['visible']]

			# Tested to thread this part as the bin_data function uses pandas functions that release the Global Interpreter Lock (GIL)
			# However, threading would only help in engineered cases where the data in each plot is huge (meaning high binning parameter but
			# even higher number of points)
			# It was decided, that the added complexity was not worth the benefit for these extreme cases. Additionally, the overhead of threads
			# would reduce performance for the majority of use cases
			binwidth = (xmax - xmin) / bins
			if len(dataframe) > max(bins, nobinning) and binwidth:
				dataframe = bin_data(dataframe, binwidth, (xmin, xmax))

			xs = dataframe['x'].to_numpy()
			ys = dataframe['y'].to_numpy() if datatype != 'lin' else np.zeros(xs.shape)

			if datatype == 'exp':
				if scaling == 'Per Plot' and len(ys):
					yrange_exp = [ys.min(), ys.max()]

				segs, colors = [], []

				filenames = dataframe['filename']
				unique_filenames = filenames.unique()
				for unique_filename in unique_filenames:
					mask = (filenames == unique_filename)
					tmp_xs, tmp_ys = xs[mask], ys[mask]

					colors.append( dataframe.loc[mask, 'color'].values )
					if cls.ids[unique_filename].is_stickspectrum:
						segs.append( np.array(((tmp_xs, tmp_xs), (np.zeros(tmp_ys.shape), tmp_ys))).T )
					else:
						segs.append( np.array(((tmp_xs[:-1], tmp_xs[1:]), (tmp_ys[:-1], tmp_ys[1:]))).T )

				if segs:
					segs = np.concatenate(segs)
					colors = np.concatenate(colors)

				self.exp_coll.set(segments=segs, colors=colors)

			elif datatype == 'cat':
				if scaling == 'Per Plot':
					yrange_cat = [ys.min(), ys.max()] if len(ys) else [-1, 1]
					ys = ys * yrange_exp[1] / yrange_cat[1]
				elif scaling in ['Global', 'Custom']:
					ys = ys * config['plot_expcat_factor'] * 10 ** config['plot_expcat_exponent']

				if config['convolution_kernel']:
					kernel_ys = np.array(config['convolution_kernel'])
					convolution_padding = len(kernel_ys) // 2
					stepwidth = config['convolution_stepwidth']

					n_xs = (xmax - xmin) // stepwidth + 1
					con_xs = np.arange(-n_xs, n_xs+1) * stepwidth + (xmax + xmin) / 2

					segs, colors = [], []

					filenames = dataframe['filename']
					unique_filenames = filenames.unique()
					for unique_filename in unique_filenames:
						mask = (filenames == unique_filename)
						tmp_xs, tmp_ys = xs[mask], ys[mask]

						hist_ys, _ = np.histogram(tmp_xs, bins=con_xs, weights=tmp_ys)
						con_ys = np.convolve(hist_ys, kernel_ys, 'full')

						con_ys = con_ys[convolution_padding:-convolution_padding]
						con_xs = con_xs[:-1]
						
						colors.append( dataframe.loc[mask, 'color'].values )
						segs.append( np.array(((con_xs[:-1], con_xs[1:]), (con_ys[:-1], con_ys[1:]))).T )

					if segs:
						segs = np.concatenate(segs)
						colors = np.concatenate(colors)

				else:
					segs = np.array(((xs, xs), (ys*0, ys))).T
					colors = dataframe['color'].to_numpy()

					mask = (xs == self.ref_position)
					colors[mask] = config['color_ref']

				self.cat_coll.set(segments=segs, colors=colors)
			
			elif datatype == 'lin':
				tuples = list(zip(xs,ys))
				tuples = tuples if len(tuples)!=0 else [[None,None]]
				colors = dataframe['color'].to_numpy()

				self.lin_coll.set_offsets(tuples)
				self.lin_coll.set_color(colors)

		if scaling == 'Per Plot':
			yrange = yrange_exp
		elif scaling == 'Global':
			yrange = ExpFile.yrange
		else:
			yrange = (config['plot_yscale_min'], config['plot_yscale_max'])

		margin = config['plot_ymargin']

		yrange = [yrange[0]-margin*(yrange[1]-yrange[0]), yrange[1]+margin*(yrange[1]-yrange[0])]
		if np.isnan(yrange[0]) or np.isnan(yrange[1]) or yrange[0] == yrange[1]:
			yrange = [-1,+1]
		ax.set_ylim(yrange)
		self.update_annotation()

	def update_annotation(self):
		fstring = config['plot_annotationfstring']

		if not fstring:
			if self.annotation:
				self.annotation.remove()
				self.annotation.set_visible(False)
				self.annotation = None
			return

		color = matplotlib.rcParams['text.color']

		if self.qns is not None:
			qns_dict = self.create_qns_dict()
			query = ' and '.join([f'({key} == {value})' for key, value in qns_dict.items()])
			
			if query and LinFile.has_results(query):
				color = config["color_lin"]
		
			noq = config['series_qns']
			qnus = [f'qnu{i+1}' for i in range(noq)]
			qnls = [f'qnl{i+1}' for i in range(noq)]

			qnus_string = ','.join([f'{qns_dict[qn]}' for qn in qnus if qn in qns_dict])
			qnls_string = ','.join([f'{qns_dict[qn]}' for qn in qnls if qn in qns_dict])

			qnstring = f'{qnus_string} ← {qnls_string}'
		else:
			qnstring = ''

		width = self.xrange[1] - self.xrange[0]

		vars = {
			'x': self.ref_position,
			'qns': qnstring,
			'width': width,
		}
		text = fstring.format(**vars)

		if self.annotation is None:
			kwargs = config['plot_annotationkwargs']
			ax = self.ax
			self.annotation = ax.text(**kwargs, s=text, color=color, transform=ax.transAxes)
		else:
			self.annotation.set_text(text)
			self.annotation.set_color(color)

	def set_xticklabels(self):
		if self.row_i:
			return
		
		ax = self.ax
		ticks = np.linspace(*self.xrange, config['plot_xticks'])
		tickformat = config['plot_xtickformat']

		if tickformat == 'absolute':
			ticklabels = symmetric_ticklabels(ticks)		
		elif tickformat == 'scientific':
			ticklabels = [f"{x:.2e}".replace("e+00", "").rstrip("0").rstrip(".") for x in ticks]
		else: # 'offset' tickformat
			ticklabels = symmetric_ticklabels(ticks - self.ref_position)

		if self.col_i and len(ticks) > 1:
			ticks = ticks[1:]
			ticklabels = ticklabels[1:]
		ax.set_xticks(ticks)
		ax.set_xticklabels(ticklabels)

	def on_range(self, xmin, xmax):
		xmin_ax, xmax_ax = self.xrange
		if xmax == xmin or xmax > xmax_ax or xmin < xmin_ax:
			return

		is_shift_pressed = (QApplication.keyboardModifiers() == Qt.KeyboardModifier.ShiftModifier)
		if is_shift_pressed: # Zoom xrange
			center, width = (xmin + xmax) / 2, xmax - xmin
			offset = center - self.ref_position

			config['plot_offset'] = offset
			config['plot_width'] = width
		else: # Fit xrange
			self.fit_data(xmin, xmax)
	
	def create_qns_dict(self, complete=False):
		qns_dict = {} if not complete else {f'qn{ul}{i+1}': pyckett.SENTINEL for ul in 'ul' for i in range(N_QNS)}
		if self.qns is None:
			return(qns_dict)
		
		qnus, qnls = self.qns
		for i, (qnu, qnl) in enumerate(zip(qnus, qnls)):
			qns_dict[f'qnu{i+1}'] = qnu
			qns_dict[f'qnl{i+1}'] = qnl
		return(qns_dict)

	@classmethod
	def fit_determine_uncert(cls, ref_pos, xmiddle, xuncert):
		error_param = config['fit_uncertainty']
		if error_param > 0:
			return(error_param)
		elif error_param >= -1:
			return( abs(xmiddle - ref_pos) )
		elif error_param >= -2:
			resp, rc = QInputDialog.getText(mainwindow, 'Set error', 'Error:')
			if rc:
				try:
					return(float(resp))
				except ValueError:
					notify_error.emit("Did not understand the given value for the error. The line was not assigned.")
					raise CustomError("Did not understand the given value for the error. The line was not assigned.")
			else:
				GUIAbortedError
		
		elif error_value >= -3:
			return(xuncert)

	def fit_data(self, xmin, xmax):
		# Delete artists highlighting previous fit
		if self.__class__.fit_vline is not None:
			self.__class__.fit_vline.remove()
			self.__class__.fit_vline = None
		if self.__class__.fit_curve is not None:
			self.__class__.fit_curve.remove()
			self.__class__.fit_curve = None
		
		# Fit the data
		xmiddle, xuncert, fit_xs, fit_ys = self.fit_peak(xmin, xmax)
		if config['fit_copytoclipboard']:
			QApplication.clipboard().setText(str(xmiddle))

		# Highlight fit in plot
		self.__class__.fit_curve = self.ax.plot(fit_xs, fit_ys, color=config["color_fit"], alpha=0.7, linewidth=1)[0]
		self.__class__.fit_vline = self.ax.axvline(x=xmiddle, color=config["color_fit"], ls="--", alpha=1, linewidth=1)

		# Create assignment object
		new_assignment = {'x': xmiddle, 'error': self.fit_determine_uncert(self.ref_position, xmiddle, xuncert), 'xpre': self.ref_position}
		new_assignment.update(self.create_qns_dict(complete=True))
		new_assignment.update({'weight': 1, 'comment': config['fit_comment'], 'filename': '__newassignments__'})
		
		if self.check_blends(new_assignment):
			return()
		NewAssignments.get_instance().add_row(new_assignment)

	def fit_peak(self, xmin, xmax):
		indices_exp = self.indices.get('exp')
		if not indices_exp:
			raise GUIAbortedError('No Experimental Indices available.')
		
		df = ExpFile.df.iloc[indices_exp[0]:indices_exp[1]]
		df = df.query(f'(visible) and x < @xmax and x > @xmin').copy()

		exp_xs, exp_ys = df['x'].to_numpy(), df['y'].to_numpy()
		fit_xs = np.linspace(xmin, xmax, config['fit_xpoints'])

		peakdirection = config['fit_peakdirection']
		fitmethod = config['fit_fitmethod']

		if (len(exp_xs) == 0) or ((len(exp_xs) < 2) and fitmethod != 'Pgopher'):
			notify_error.emit('The data could not be fit as there were too few points selected.')
			raise GUIAbortedError('The data could not be fit as there were too few points selected.')

		try:
			fit_function = get_fitfunction(fitmethod, config['fit_offset'])
			xmiddle, xuncert, fit_xs, fit_ys = fit_function(exp_xs, exp_ys, peakdirection, fit_xs)
		except Exception as E:
			self.fitcurve = None
			self.fitline = None
			notify_error.emit(f"The fitting failed with the following error message : {str(E)}")
			raise

		return(xmiddle, xuncert, fit_xs, fit_ys)

	def check_blends(self, new_assignment):
		if not config['series_blendwidth']:
			return(False)
		reference_states = ReferenceSeriesWindow.instance.get_state()
		if self.col_i > len(reference_states):
			return(False)
		reference_state = reference_states[self.col_i]
		if not reference_state['check_blends']:
			return(False)

		
		if reference_state['method'] == 'Transition':
			file_to_limit_data_to = reference_state['transition'].get('file')
		else:
			file_to_limit_data_to = False
		
		blendwidth = config['series_blendwidth']
		xpre = new_assignment['xpre']
		xmin, xmax = xpre - blendwidth, xpre + blendwidth
		
		
		entries = CatFile.get_data(xrange=(xmin, xmax)).copy()
		if file_to_limit_data_to:
			query = 'filename == @file_to_limit_data_to'
			entries = entries.query(query)
		
		if config['series_blendminrelratio']:
			qn_labels = pyckett.qnlabels_from_quanta(N_QNS)
			query = ' and '.join([f'( {label} == {new_assignment[label]} )' for label in qn_labels if new_assignment[label] != pyckett.SENTINEL])
			tmp_cat = entries.query(query)
			y_at_xpre = tmp_cat["y"].values[0]
		
			min_y_value = y_at_xpre * config['series_blendminrelratio']
			entries = entries.query('y >= @min_y_value')
			
		if len(entries) < 2:
			return(False)
		
		# After guard clauses we want to show the dialog
		mainwindow.lwpwidget.drawplot.emit()
		AssignBlendsDialog.show_dialog(new_assignment, entries)
		return(True)

class LWPWidget(QGroupBox):
	drawplot = pyqtSignal()
	plotscreated = pyqtSignal()
	draw_counter = AtomicCounter()
	_active_ax_index = (0, 0)

	_ax_class = LWPAx

	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		self.initialize_dynamic_decorator_for_drawing_plot()

		with matplotlib_lock:
			self.fig = matplotlib.figure.Figure(dpi=config["plot_dpi"])
			self.listener_onclick = self.fig.canvas.mpl_connect('button_press_event', self.on_click)
			self.listener_onhover = self.fig.canvas.mpl_connect("motion_notify_event", self.on_hover)
			self.plotcanvas = FigureCanvas(self.fig)
			self.plotcanvas.contextMenuEvent = self.contextMenuCanvas
	
		config.register("plot_dpi", lambda: self.fig.set_dpi(config["plot_dpi"]))
		
		self.plotcanvas.setMinimumHeight(200)
		self.plotcanvas.setMinimumWidth(200)
		
		self.drawplot.connect(self.draw_canvas)
		keys = ('plot_width', 'plot_widthexpression', 'plot_offset', 'plot_offsetexpression', 'plot_offsetisrelative', 'plot_bins',
				'plot_yscale', 'plot_yscale_min', 'plot_yscale_max', 'plot_expcat_factor', 'plot_expcat_exponent')
		config.register(keys, lambda: self.set_data())

		self.mpltoolbar = NavigationToolbar2QT(self.plotcanvas, self)
		self.mpltoolbar.setVisible(config["isvisible_matplotlibtoolbar"])
		config.register('isvisible_matplotlibtoolbar', lambda: self.mpltoolbar.setVisible(config['isvisible_matplotlibtoolbar']))

		toplayout = self.toplayout = QHBoxLayout()

		buttonsdict = {
			'in':			lambda _: WidthDialog.gui_set_width('++'),
			'out':			lambda _: WidthDialog.gui_set_width('--'),
			'left':			lambda _: OffsetDialog.gui_set_offset('-'),
			'right':		lambda _: OffsetDialog.gui_set_offset('+'),
		}

		for label, func in buttonsdict.items():
			button = QQ(QPushButton, text=label, change=func, visible=config['isvisible_controlsmainplot'])
			toplayout.addWidget(button)
			config.register('isvisible_controlsmainplot', lambda button=button: button.setVisible(config['isvisible_controlsmainplot']))

		toplayout.addStretch(1)

		rows_cols_elements = (
			QQ(QLabel, text="Plots: ", visible=config["isvisible_controlsrowscols"]),
			QQ(QSpinBox, "plot_rows", range=(1, None), maxWidth=45, visible=config["isvisible_controlsrowscols"]),
			QQ(QLabel, text="x", visible=config["isvisible_controlsrowscols"]),
			QQ(QSpinBox, "plot_cols", range=(1, None), maxWidth=45, visible=config["isvisible_controlsrowscols"]),
		)
		for elem in rows_cols_elements:
			toplayout.addWidget(elem)
			config.register("isvisible_controlsrowscols", lambda elem=elem: elem.setVisible(config["isvisible_controlsrowscols"]))

		width_elements = (
			QQ(QLabel, text="    Width: ", visible=config["isvisible_controlswidth"]),
			QQ(QDoubleSpinBox, "plot_width", range=(0, None), minWidth=85, visible=config["isvisible_controlswidth"]),
		)
		for elem in width_elements:
			toplayout.addWidget(elem)
			config.register("isvisible_controlswidth", lambda elem=elem: elem.setVisible(config["isvisible_controlswidth"]))

		layout = QVBoxLayout()
		layout.addLayout(toplayout)
		layout.addWidget(self.plotcanvas, 1)
		layout.addWidget(self.mpltoolbar)
		self.setLayout(layout)

		self.lwpaxes = np.zeros((0, config['plot_cols']))

		config.register(("plot_rows", "plot_cols"), self.create_plots)
		self.create_plots()

	def initialize_dynamic_decorator_for_drawing_plot(self):
		def exit_func():
			counter_value = self.draw_counter.decrease()
			if not counter_value:
				self.drawplot.emit()

		def init_func():
			self.draw_counter.increase()

		drawplot_decorator.init_func = init_func
		drawplot_decorator.exit_func = exit_func

	def draw_canvas(self):
		with matplotlib_lock:
			self.plotcanvas.draw_idle()

	def on_click(self, event):
		ax = event.inaxes
		index = np.asarray(np.where(self.axes  == ax)).T
		
		if len(index):
			self.__class__._active_ax_index = tuple(index[0])

	def on_hover(self, event):
		x = event.xdata
		y = event.ydata

		if not all([x, y, event.inaxes]):
			text_cursor = ""
		else:
			if config['flag_showmainplotposition']:
				text_cursor = f"  ({x=:{config['flag_xformatfloat']}}, {y=:{config['flag_xformatfloat']}})  "
			else:
				text_cursor = ""

			CloseByLinesWindow.instance.cursor_changed.emit(x, y)
		mainwindow.statusbar.position_label.setText(text_cursor)

	def wheelEvent(self,event):
		steps = event.angleDelta().y() // 120
		factor = 2**(-steps)
		WidthDialog.gui_set_width(factor)

	@QThread.threaded_d
	@lock_d(matplotlib_lock)
	@status_d
	def create_plots(self, thread=None):
		n_rows = max(config["plot_rows"], 1)
		n_cols = max(config["plot_cols"], 1)
		gridspec_kw = config["plot_gridspeckwargs"]

		# The following performed better than self.fig.clf()
		for ax in self.fig.get_axes():
			self.fig.delaxes(ax)
		
		for lwpax in self.lwpaxes.flatten():
			lwpax.span = None

		thread.earlyreturn()

		axes = self.fig.subplots(n_rows, n_cols, gridspec_kw=gridspec_kw, squeeze=False)[::-1]
		self.axes = axes

		thread.earlyreturn()

		indices = np.indices(axes.shape)
		vectorized_ax_class = np.vectorize(self._ax_class, otypes=[object])
		self.lwpaxes = vectorized_ax_class(axes, *indices)

		thread.earlyreturn()

		self.set_data()
		self.plotscreated.emit()
	
	def calculate_widths(self):
		expression = config['plot_widthexpression']
		value = config['plot_width']
		if not expression:
			return(value)
		
		compiled_expression = compile(expression, '', 'eval')
		return(compiled_expression)
		
	@QThread.threaded_d
	@status_d
	@drawplot_decorator.d
	def set_data(self, thread=None):
		if not hasattr(ReferenceSeriesWindow, 'instance'):
			return
		
		if not hasattr(ReferenceSeriesWindow.instance, 'tab'):
			return

		n_rows = config['plot_rows']
		n_cols = config['plot_cols']
		n_qns = config['series_qns']

		# Setup arrays to hold final data
		positions_shape = (n_rows, n_cols)
		positions = np.empty(positions_shape, dtype=np.float64)
		qns_shape = (n_rows, n_cols, 2, n_qns)
		qns = np.empty(qns_shape, dtype=np.int64)

		thread.earlyreturn()

		# Calculate positions and qns for each column
		tab_widget = ReferenceSeriesWindow.instance.tab
		n_widgets = tab_widget.count()

		for i_col in range(n_widgets):
			thread.earlyreturn()
			if i_col > n_cols:
				continue

			refwidget = tab_widget.widget(i_col)
			positions[:, i_col], qns[:, i_col] = refwidget.calc_references(n_rows, n_qns)
		
		thread.earlyreturn()

		# Edge cases of no reference tabs and too few tabs
		if n_widgets < 1:
			positions[:] = 0
			qns[:] = 0

		elif n_widgets < n_cols:
			for i_col in range(n_widgets, n_cols):
				positions[:, i_col] = positions[:, n_widgets - 1]
				qns[:, i_col] = qns[:, n_widgets - 1]

		thread.earlyreturn()

		# Calculate the widths and offsets
		width_value = config['plot_width']
		width_expr = config['plot_widthexpression']

		offset_value = config['plot_offset']
		offset_expr = config['plot_offsetexpression']

		if not width_expr:
			widths = width_value
		else:
			width_expr_comp = compile(width_expr, '', 'eval')
			i_rows, i_cols = np.indices((n_rows, n_cols))
			vars = {'x': positions, 'w': width_value, 'i_row': i_rows, 'i_col': i_cols}
			for i in range(n_qns):
				vars[f'qnu{i+1}'] = qns[:,:,0,i]
				vars[f'qnl{i+1}'] = qns[:,:,1,i]

			widths = eval(width_expr_comp, vars)
		
		if not offset_expr:
			offsets = offset_value
		else:
			offset_expr_comp = compile(offset_expr, '', 'eval')
			i_rows, i_cols = np.indices((n_rows, n_cols))
			vars = {'x': positions, 'o': offset_value, 'i_row': i_rows, 'i_col': i_cols}
			for i in range(n_qns):
				vars[f'qnu{i+1}'] = qns[:,:,0,i]
				vars[f'qnl{i+1}'] = qns[:,:,1,i]

			offsets = eval(offset_expr_comp, vars)

		thread.earlyreturn()

		# Calculate the plot positions and get indices
		xmins = positions + offsets - widths/2
		xmaxs = xmins + widths

		min_indices, max_indices = {}, {}
		for label, cls in {'exp': ExpFile, 'cat': CatFile, 'lin': LinFile}.items():
			min_indices[label], max_indices[label] = cls.xs_to_indices(xmins, xmaxs)
		
		thread.earlyreturn()

		# Set the correct values to the LWPAxes
		if self.lwpaxes.shape != (n_rows, n_cols):
			notify_error.emit('Shape of LWPAxes is out of sync with requested values.')
			return
		
		threads = []
		for i_row in range(n_rows):
			for i_col in range(n_cols):
				ax = self.lwpaxes[i_row, i_col]
				ax.ref_position = positions[i_row, i_col]
				ax.xrange = (xmins[i_row, i_col], xmaxs[i_row, i_col])
				ax.qns = qns[i_row, i_col]
				ax.indices = {label: (min_indices[label][i_row, i_col], max_indices[label][i_row, i_col]) for label in ('exp', 'cat', 'lin')}
				threads.append(ax.update())

		thread.earlyreturn()

		for thread_ in threads:
			thread_.wait()

	def get_current_ax(self):
		shape = self.lwpaxes.shape
		index = self.__class__._active_ax_index
		if 0 <= index[0] < shape[0] and 0 <= index[1] < shape[1]:
			return(self.lwpaxes[index])
		else:
			return(self.lwpaxes[0,0])

	def contextMenuCanvas(self, event):
		x, y = event.x(), event.y()
		geometry = self.plotcanvas.geometry()
		width, height = geometry.width(), geometry.height()
		x_rel, y_rel = x/width, 1 - y/height

		for lwpax in self.lwpaxes.flatten():
			xmin, ymin, width, height = lwpax.ax.get_position().bounds
			if xmin <= x_rel <= xmin + width and ymin <= y_rel <= ymin+height:
				break
		else: # Clicked outside of ax
			return

		menu = QMenu(self)
		get_position_action = menu.addAction('Copy Reference Position')
		get_qns_action = menu.addAction('Copy QNs')
		set_active_action = menu.addAction('Make active Ax')
		fit_all_action = menu.addAction('Fit all')

		action = menu.exec(self.mapToGlobal(event.pos()))
		if action == get_position_action:
			QApplication.clipboard().setText(str(lwpax.ref_position))
		elif action == get_qns_action:
			output_string = []
			qnus, qnls = lwpax.qns
			tmp_string = ','.join([f'{qn:3.0f}' for qn in qnus]) + ' ← ' + ','.join([f'{qn:3.0f}' for qn in qnls])
			output_string.append(tmp_string)

			output_string = '\n'.join(output_string)
			QApplication.clipboard().setText(output_string)
		elif action == set_active_action:
			mainwindow.lwpwidget._active_ax_index = (lwpax.row_i, lwpax.col_i)
		elif action == fit_all_action:
			i_col = lwpax.col_i
			AssignAllDialog.show_dialog(i_col)

class Menu():
	def __init__(self, parent, *args, **kwargs):
		mb = parent.menuBar()
		
		# Create top level menus
		top_menu_labels = ("Files", "View", "Fit", "Plot", "Modules", "Info")
		self.top_menus = {}
		for label in top_menu_labels:
			menu = mb.addMenu(f'{label}')
			self.top_menus[label] = menu
		

		toggleaction_files = FileWindow.instance.toggleViewAction()
		toggleaction_files.setText('Edit Files')
		toggleaction_files.setShortcut('Shift+1')

		toggleaction_config = ConfigWindow.instance.toggleViewAction()
		toggleaction_config.setShortcut('Shift+0')

		toggleaction_convolution = ConvolutionWindow.instance.toggleViewAction()
		toggleaction_convolution.setText('Sticks to Lineshape')
		toggleaction_convolution.setToolTip('Choose a function to create a spectrum from stick data')

		toggleaction_credits = CreditsWindow.instance.toggleViewAction()
		toggleaction_credits.setText("Credits and License")
		toggleaction_credits.setToolTip("See the Credits and License")

		view_actions = [		
			ReferenceSeriesWindow.instance.toggleViewAction(),
			NewAssignmentsWindow.instance.toggleViewAction(),
			CloseByLinesWindow.instance.toggleViewAction(),
			LogWindow.instance.toggleViewAction(),
		]
		
		for i, view_action in enumerate(view_actions):
			view_action.setShortcut(f'Shift+{i+2}')

		modules_actions = [
			ResidualsWindow.instance.toggleViewAction(),
			BlendedLinesWindow.instance.toggleViewAction(),
			ReportWindow.instance.toggleViewAction(),
			SeriesfinderWindow.instance.toggleViewAction(),
			PeakfinderWindow.instance.toggleViewAction(),
			EnergyLevelsWindow.instance.toggleViewAction(),
			CmdWindow.instance.toggleViewAction(),
		]

		for i, modules_action in enumerate(modules_actions):
			modules_action.setShortcut(f'Ctrl+{i+1}')

		fitfunction_menu = QMenu("Choose Fit Function", parent=parent)
		self.fitfunction_actions = {}

		current_method = config['fit_fitmethod']
		for method in LWPAx.fit_methods:
			is_checked = (method == current_method)
			callback = lambda _, method=method: self.set_fitmethod(method)
			self.fitfunction_actions[method] = QQ(QAction, parent=parent, text=f"{method}", change=callback, checkable=True, value=is_checked)
			fitfunction_menu.addAction(self.fitfunction_actions[method])
		config.register('fit_fitmethod', self.on_fitfunction_changed)

		actions = {
			'Files': (
				QQ(QAction, parent=parent, text="Add Files", change=File.add_files_dialog, shortcut='Ctrl+O', tooltip="Add any kind of Files"),
				QQ(QAction, parent=parent, text='Reread All Files', change=lambda _: File.reread_all(), shortcut='Ctrl+R', tooltip="Reread all Exp, Cat and Lin files"),
				QQ(QAction, 'flag_autoreloadfiles', checkable=True, parent=parent, text='Auto Reload Files', tooltip="Automatically reload files on change"),
				None,
				toggleaction_files,
				None,
				QQ(QAction, parent=parent, text="Save current values as default", tooltip="Save current configuration as default", change=lambda _: config.save()),
				None,
				QQ(QAction, parent=parent, text="Save Files as Project", change=lambda _: File.save_files_gui(), tooltip="Save all loaded files and their parameters as a project."),
				QQ(QAction, parent=parent, text="Load Project", change=lambda _: File.load_files_gui(), tooltip="Load a project."),
				None,
				QQ(QAction, parent=parent, text="Quit", change=mainwindow.close),
			),
			'Fit': (
				fitfunction_menu,
				QQ(QAction, parent=parent, text="Change Function", shortcut="Ctrl+F", tooltip="Cycle through the available fit-functions", change=lambda _: self.next_fitmethod()),
				None,
				QQ(QAction, parent=parent, text="Change Fit Color", tooltip="Change the color of the fitfunction", change=lambda _: self.change_fitcolor()),
			),
			'Plot': (
				QQ(QAction, parent=parent, text="# Plots", shortcut="Ctrl+N", tooltip="Change number of plots", change=lambda _: self.plot_number()),
				None,
				QQ(QAction, parent=parent, text="Set Width", shortcut="Ctrl+W", change=lambda _: WidthDialog.show_dialog()),
				QQ(QAction, parent=parent, text="Set Offset", shortcut="Ctrl+G", change=lambda _: OffsetDialog.show_dialog()),
				None,
				ScalingWindow.instance.toggleViewAction(),
				toggleaction_convolution,
			),
			'View': (
				toggleaction_config,
				None,
				*view_actions,
				None,
				QQ(QAction, parent=parent, text='Open Console', shortcut='CTRL+K', change= lambda _: ConsoleDialog.show_dialog()),
				None,
			),
			'Modules': modules_actions,
			'Info': (
				QQ(QAction, parent=parent, text="Open LLWP folder", change=lambda x: webbrowser.open(f'file:///{llwpfile()}'), tooltip="Open the folder containing the config, ...", ),
				QQ(QAction, parent=parent, text="Send Mail to Author", change=lambda x: self.send_mail_to_author(), tooltip="Send a mail to the developer"),
				toggleaction_credits,
			)
			
		}

		for label, menu in self.top_menus.items():
			for widget in actions.get(label, []):
				if widget is None:
					menu.addSeparator()
				elif isinstance(widget, QAction):
					menu.addAction(widget)
				else:
					menu.addMenu(widget)
		
	def next_fitmethod(self):
		fitmethods = LWPAx.fit_methods
		newindex = fitmethods.index( config['fit_fitmethod'] ) + 1
		newindex = newindex % len(fitmethods)
		newvalue = fitmethods[newindex]
		config['fit_fitmethod'] = newvalue

	def set_fitmethod(self, method):
		config['fit_fitmethod'] = method

	def on_fitfunction_changed(self, nextfunction=False):
		value = config["fit_fitmethod"]
		for fitfunction, action in self.fitfunction_actions.items():
			action.setChecked(fitfunction == value)
		notify_info.emit(f"Fitting with method '{value}'")

	def change_fitcolor(self):
		color = QColorDialog.getColor(initial=QColor(Color.rgbt_to_trgb(config['color_fit'])), options=QColorDialog.ColorDialogOption.ShowAlphaChannel)
		color = Color(Color.trgb_to_rgbt(color.name(QColor.NameFormat.HexArgb)))
		config['color_fit'] = color
	
	def plot_number(self):
		resp, rc = QInputDialog.getText(mainwindow, 'How many plots do you want: ', 'Number:')
		if not rc:
			return
		resp = resp.split("x")
		try:
			if len(resp) == 1:
				config["plot_rows"] = int(resp[0])
			elif len(resp) == 2:
				config["plot_rows"] = int(resp[0])
				config["plot_cols"] = int(resp[1])
		except ValueError:
			notify_warning.emit("The entered value was not understood.")
			return
		
	def send_mail_to_author(self):
		webbrowser.open(f"mailto:bonah@ph1.uni-koeln.de?subject={APP_TAG}")


class MainWindow(QMainWindow):
	notify_info = pyqtSignal(str)
	notify_warning = pyqtSignal(str)
	notify_error = pyqtSignal(str)

	status_counter = AtomicCounter()
	working = pyqtSignal()
	waiting = pyqtSignal()

	mainwidget_class = LWPWidget
	menu_class = Menu

	def __init__(self, app, *args, **kwargs):
		super().__init__(*args, **kwargs)
		self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
		self.setWindowTitle(APP_TAG)
		self.setAcceptDrops(True)
		Geometry.load_widget_geometry(self)
		
		try:
			# Set LLWP logo as icon
			possible_folders = [os.path.dirname(os.path.realpath(__file__)), os.getcwd()]
			for folder in possible_folders:
				iconpath = os.path.join(folder, "LLWP.svg")
				if os.path.isfile(iconpath):
					icon = QIcon(iconpath)
					break
			
			# Make LLWP appear separate from python scripts in taskbar
			app.setWindowIcon(QIcon(iconpath))
			import ctypes
			ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(APP_TAG)
		except Exception as E:
			pass
		
	def create_gui_components(self):
		global notify_info, notify_warning, notify_error
		notify_info = self.notify_info
		notify_warning = self.notify_warning
		notify_error = self.notify_error


		self.statusbar = StatusBar()
		self.setStatusBar(self.statusbar)

		self.notificationsbox = NotificationsBox()

		mainwidget_class = self.__class__.mainwidget_class
		self.lwpwidget = mainwidget_class()
		self.setCentralWidget(self.lwpwidget)
				
		NewAssignments.get_instance()
		
		for subclass in get_all_subclasses(EQDockWidget):
			if APP_TAG in subclass.available_in:
				subclass()
		dockstate = Geometry.get('__dockstate__')
		if dockstate:
			self.restoreState(dockstate)
		
		menu_class = self.__class__.menu_class
		self.menu = menu_class(self)
		self.shortcuts()

	def shortcuts(self):
		shortcuts_dict = {
			"w": lambda: WidthDialog.gui_set_width("++"),
			"s": lambda: WidthDialog.gui_set_width("--"),
			"a": lambda: OffsetDialog.gui_set_offset("--"),
			"d": lambda: OffsetDialog.gui_set_offset("++"),

			"Shift+w": lambda: WidthDialog.gui_set_width("+"),
			"Shift+s": lambda: WidthDialog.gui_set_width("-"),
			"Shift+a": lambda: OffsetDialog.gui_set_offset("-"),
			"Shift+d": lambda: OffsetDialog.gui_set_offset("+"),
			
			"Ctrl+S": lambda: NewAssignments.get_instance().save_gui(),
			"Ctrl+Space": lambda: AssignAllDialog.show_dialog(),
			"Ctrl+Shift+k": lambda: ConsoleDialog.run_current_command(),
			"F11": lambda: self.togglefullscreen(),
		}

		for keys, function in shortcuts_dict.items():
			tmp = QShortcut(keys, self)
			tmp.activated.connect(function)
			tmp.setContext(Qt.ShortcutContext.WidgetShortcut)

	def togglefullscreen(self):
		if self.isFullScreen():
			self.showNormal()
		else:
			self.showFullScreen()
	

	def dragEnterEvent(self, event):
		if event.mimeData().hasUrls():
			event.accept()
		else:
			event.ignore()

	def dropEvent(self, event):
		event.accept()
		files_by_type = File.sort_files_by_type([url.toLocalFile() for url in event.mimeData().urls()])
		File.add_multiple_files_by_type(files_by_type)

	def moveEvent(self, *args, **kwargs):
		Geometry.save_widget_geometry(self)
		return super().moveEvent(*args, **kwargs)

	def resizeEvent(self, *args, **kwargs):
		Geometry.save_widget_geometry(self)
		return super().resizeEvent(*args, **kwargs)

	def closeEvent(self, *args, **kwargs):
		Geometry.set("__dockstate__", self.saveState())
		Geometry.save()
		config.save()


class StatusBar(QStatusBar):
	set_cursor_text = pyqtSignal(str)

	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)

		self._disappearing_messages = []
		self._disappearing_messages_timer = QTimer()
		self._disappearing_messages_timer.setSingleShot(True)
		self._disappearing_messages_timer.timeout.connect(self.next_disappearing_message)

		self.setStyleSheet("QStatusBar::item {border-left: 1px solid inherit;}")

		self.messages_label = QQ(QLabel, text='', wordwrap=True)
		self.addPermanentWidget(self.messages_label, 1)

		self.position_label = QQ(QLabel, text='')
		self.addPermanentWidget(self.position_label)

		self.working_label = QQ(QLabel, text='')
		self.addPermanentWidget(self.working_label)
		
		mainwindow.working.connect(lambda: self.working_label.setText("   Working ...  "))
		mainwindow.waiting.connect(lambda: self.working_label.setText("   Ready  "))

		notify_error.connect(lambda x: self.disappearing_message(x, style='error'))
		notify_warning.connect(lambda x: self.disappearing_message(x, style='warning'))
		notify_info.connect(lambda x: self.disappearing_message(x, style='info'))

	def disappearing_message(self, text, style=None):
		style_string = {
			'error': "<span style='color:#ff0000;'>ERROR</span>: ",
			'warning': "<span style='color:#eda711;'>WARNING</span>: ",
			'info': "<span style='color:#0096FF;'>INFO</span>: ",
		}.get(style, '')

		max_characters = config['flag_statusbarmaxcharacters']
		
		if len(text) > max_characters:
			text = f'{text[:max_characters]} ...'
		
		text = f'  {style_string}{text}'
		self._disappearing_messages.append(text)

		if not self._disappearing_messages_timer.isActive():
			self.next_disappearing_message()

	def next_disappearing_message(self):
		try:
			if len(self._disappearing_messages) > 3:
				self._disappearing_messages = []
				next_message = 'Many new messages. Please see the log window for all messages.'
			else:
				next_message = self._disappearing_messages.pop(0)
		except IndexError:
			self.messages_label.setText('')
			return
		
		self.messages_label.setText(next_message)
		self._disappearing_messages_timer.start(config['flag_notificationtime'])

class NotificationsBox(QScrollArea):
	layout_updated = pyqtSignal()

	def __init__(self):
		super().__init__()
		
		self.setWindowFlags(
			Qt.WindowType.Window | Qt.WindowType.Tool | Qt.WindowType.FramelessWindowHint |
			Qt.WindowType.WindowStaysOnTopHint | Qt.WindowType.X11BypassWindowManagerHint)
		
		self.setAttribute(Qt.WidgetAttribute.WA_NoSystemBackground, True)
		self.setAttribute(Qt.WidgetAttribute.WA_TranslucentBackground, True)
		
		self.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)

		self.messages = []

		self.widget = QWidget()
		self.setWidget(self.widget)
		self.layout = QVBoxLayout()
		self.widget.setLayout(self.layout)

		available_geometry = QApplication.primaryScreen().availableGeometry()
		self.available_height = available_geometry.height()
		self.available_width = available_geometry.width()

		self.padding = 10
		self.width_ = 300
		window_position = QPoint(self.available_width - self.width_ - self.padding, self.padding)
		self.move(window_position)

		self.max_height = self.available_height - 2 * self.padding
		self.setFixedHeight(self.max_height)
		self.setFixedWidth(self.width_)


		self.setStyleSheet(f"""
			background-color: transparent;
		""")


		self.widget.setStyleSheet("""
			background-color: transparent;
		""")

		self.layout_updated.connect(self.update_layout)

		notify_error.connect(lambda x: self.add_message(x, style='error'))
		notify_warning.connect(lambda x: self.add_message(x, style='warning'))
		notify_info.connect(lambda x: self.add_message(x, style='info'))


	def add_message(self, text, style=None):
		style_string = {
			'error': "<span style='color:#ff0000;'>ERROR</span>: ",
			'warning': "<span style='color:#eda711;'>WARNING</span>: ",
			'info': "<span style='color:#0096FF;'>INFO</span>: ",
		}.get(style, '')

		text = f'{style_string}{text}'

		label = QQ(QLabel, text=text, width=self.width_-20, wordwrap=True)
		label.setStyleSheet("""
			padding: 5px;
			border-radius: 5px;
			color: white;
			background-color: #bf29292a;
		""")

		self.layout.addWidget(label)
		self.messages.append(label)

		self.show()

		disappear_timer = QTimer(self)
		disappear_timer.setSingleShot(True)
		disappear_timer.timeout.connect(self.unshow_message)
		disappear_timer.start(config["flag_notificationtime"])
		self.layout_updated.emit()

	def update_layout(self):
		layout_timer = QTimer(self)
		layout_timer.setSingleShot(True)
		layout_timer.timeout.connect(lambda: [self.widget.adjustSize(), self.adjustSize()])
		layout_timer.start(0)

	def unshow_message(self):
		label = self.messages.pop()
		label.hide()
		label.deleteLater()
		if not self.messages:
			self.hide()
		self.layout_updated.emit()

class LLWP(QApplication):
	configsignal = pyqtSignal(tuple)
	mainwindow_class = MainWindow
	
	def __init__(self, *args, **kwargs):
		sys.excepthook = except_hook
		threading.excepthook = lambda args: except_hook(*args[:3])
		
		super().__init__(sys.argv, *args, **kwargs)
		self.setStyle("Fusion")

		Geometry().load()

		global config
		config = Config(self.configsignal)
		config.load()
		messages = config.messages
		
		config.register('flag_pyckettquanta', self.update_pyckett_quanta_settings)
		self.update_pyckett_quanta_settings()

		self.initialize_matplotlib_settings()

		with config:
			global mainwindow
			mainwindow = self.mainwindow_class(self)
			mainwindow.create_gui_components()

			if len(sys.argv) > 1:
				project = sys.argv[1]
				File.load_files(project)

			mainwindow.show()
			
			if messages:
				notify_warning.emit('\n'.join(messages))
			
			self.run_init_commands()
			self.debug_setup()

		sys.exit(self.exec())
	
	# Used to automatically reproduce behavior when testing
	def debug_setup(self):
		pass

	def run_init_commands(self):
		try:
			commands = config['commandlinedialog_commands']
			for _, command, run_initially in commands:
				if run_initially:
					ConsoleDialog.run_command(command)
		except Exception as E:
			notify_error.emit(f'{E}:\n{tb.format_stack()}')

	def initialize_matplotlib_settings(self):
		matplotlib.rc('font', **config["plot_fontdict"])
		config.register('plot_fontdict', lambda: matplotlib.rc('font', **config['plot_fontdict']))
		self.styleHints().colorSchemeChanged.connect(self.update_matplotlib_theme)
		self.update_matplotlib_theme()

	def update_matplotlib_theme(self):
		matplotlib.style.use('dark_background' if is_dark_theme() else 'default')
		matplotlib.rcParams['figure.facecolor'] = '#00000000'
		matplotlib.rcParams['axes.facecolor'] = '#00000000'

	def update_pyckett_quanta_settings(self):
		pyckett.QUANTA = config['flag_pyckettquanta']


##
## Dialogs
##
class AssignBlendsDialog(QDialog):
	_instance = None

	@classmethod
	def show_dialog(cls, new_assignment, entries):
		dialog = cls(new_assignment, entries)
		dialog.show()

	def __init__(self, new_assignment, entries):
		super().__init__()
		self.setModal(False)
		
		self.setWindowTitle("Assign Blends")

		self.finished.connect(self.on_exit)
		if self._instance:
			self._instance.done(0)
		self.__class__._instance = self
		

		self.init_gui()
		self.update_gui(new_assignment, entries)

	def update_gui(self, new_assignment, entries):
		self.noq = config["series_qns"]
		self.entries = entries.reset_index(drop=True)
		self.xmiddle = new_assignment['x']
		self.error = new_assignment['error']

		self.label.setText(f'Assigning Blend for position {self.xmiddle}.')

		self.entries["dist"] = self.entries["x"] - self.xmiddle
		self.entries["absdist"] = abs(self.entries["dist"])
		self.entries.sort_values(by=["absdist"], inplace=True)
		self.entries["y"] = np.log(self.entries["y"])/np.log(10)

		noq = config["series_qns"]
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
					text = f"{entry[val]:{config['flag_xformatfloat']}}".rstrip("0").rstrip(".")
				table.setItem(currRowCount, j+1, QTableWidgetItem(text))
			checkbox = QQ(QCheckBox, value=True)
			self.checkboxes.append(checkbox)
			table.setCellWidget(currRowCount, 0, checkbox)

		for i in range(N_QNS):
			table.setColumnHidden(i+ 4, i>=noq)
			table.setColumnHidden(i+ 4 + N_QNS, i>=noq)

		table.resizeColumnsToContents()

	def init_gui(self):
		layout = QVBoxLayout(margin=True)
		self.setLayout(layout)

		self.label = QQ(QLabel)
		layout.addWidget(self.label)

		self.table = QTableWidget()
		self.cols = ["x", "y", "dist"] + [f"qn{ul}{i+1}" for ul in ("u", "l") for i in range(N_QNS)] + ["filename"]
		self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
		self.table.setRowCount(0)
		self.table.setColumnCount(len(self.cols)+1)
		self.table.setHorizontalHeaderLabels(["Y/N", "x", "log. y", "Dist"] +  [f"{ul}{i+1}" for ul in ("U", "L") for i in range(N_QNS)] + ["Filename"])
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
		self.done(0)
		checked = True if all else [x.isChecked() for x in self.checkboxes]
		self.entries['checked'] = checked

		selected_rows = self.entries.query('checked == True')
		max_predicted = selected_rows['y'].max()

		new_assignments = {
			'x': self.xmiddle,
			'weight': 10 ** (selected_rows['y'] - max_predicted),
			'error': self.error,
			'comment': config['fit_comment'],
			'filename': '__newassignments__',
		}

		for i in range(self.noq):
			new_assignments[f'qnu{i+1}'] = selected_rows[f'qnu{i+1}']
			new_assignments[f'qnl{i+1}'] = selected_rows[f'qnl{i+1}']
		
		for i in range(self.noq, N_QNS):
			new_assignments[f'qnu{i+1}'] = pyckett.SENTINEL
			new_assignments[f'qnl{i+1}'] = pyckett.SENTINEL

		NewAssignments.get_instance().add_rows(new_assignments)

	def on_exit(self):
		self.__class__._instance = None

class ConsoleDialog(QDialog):
	def __init__(self):
		super().__init__()
		self.setWindowTitle('Command Line Dialog')

		self.tabs = QTabWidget()
		self.tabs.setTabsClosable(True)
		self.tabs.setMovable(True)
		self.tabs.setDocumentMode(True)

		initial_values = config['commandlinedialog_commands']
		if initial_values:
			for title, command, run_initially in initial_values:
				self.add_tab(title, command, run_initially)
		else:
			self.add_tab()

		self.tabs.tabCloseRequested.connect(self.close_tab)
		self.tabs.tabBarDoubleClicked.connect(self.renameoradd_tab)
		self.tabs.setCurrentIndex(config['commandlinedialog_current'])

		layout = QVBoxLayout(margin=True)
		self.setLayout(layout)

		layout.addWidget(self.tabs)
		buttons_layout = QHBoxLayout()
		buttons_layout.addStretch()
		buttons_layout.addWidget(QQ(QPushButton, text='Run', change=lambda x: self.done(1), shortcut='Ctrl+Return'))
		buttons_layout.addWidget(QQ(QPushButton, text='Cancel', change=lambda x: self.done(0), shortcut='Esc'))
		buttons_layout.addStretch()
		layout.addLayout(buttons_layout)

	def add_tab(self, title='Command', command='', run_initially=False):
		widget = QWidget()
		layout = QVBoxLayout()
		widget.setLayout(layout)

		textarea = QQ(QPlainTextEdit, value=command)
		cursor = textarea.textCursor()
		cursor.movePosition(QTextCursor.MoveOperation.End)
		textarea.setTextCursor(cursor)

		runinit_checkbox = QQ(QCheckBox, value=run_initially, text='Run on startup')

		layout.addWidget(textarea)
		layout.addWidget(runinit_checkbox)
		self.tabs.addTab(widget, title)

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
			text, ok = QInputDialog().getText(self, "Tab Name", "Enter the Tabs Name:")
			if ok and text:
				self.tabs.setTabText(index, text)

	def done(self, val):
		super().done(val)
		commands = []
		for i in range(self.tabs.count()):
			tab = self.tabs.widget(i)
			title = self.tabs.tabText(i)

			command = tab.findChildren(QPlainTextEdit)[0].toPlainText()
			run_initially = tab.findChildren(QCheckBox)[0].isChecked()
			commands.append((title, command, run_initially))

		config['commandlinedialog_commands'] = commands
		config['commandlinedialog_current'] = self.tabs.currentIndex()

	
	@staticmethod
	def run_current_command():
		commands = config['commandlinedialog_commands']
		if not commands:
			notify_warning.emit('No commands specified. Therfore no command was executed.')
			return
	
		index = config['commandlinedialog_current']
		if index >= len(commands):
			notify_warning.emit('The index of the command is out of the range.')
			return
		
		_, command, _ = commands[index]
		if not command.strip():
			return
		ConsoleDialog.run_command(command)
		
	
	@staticmethod
	def run_command(command):
		message = []
		old_stdout = sys.stdout
		red_output = sys.stdout = io.StringIO()
		try:
			exec(command)
		except Exception as E:
			message.append(f'Executing the code raised an error: {str(E)}')
			raise
		finally:
			sys.stdout = old_stdout

		message.append('\n'.join([f'>>> {line}' for line in command.split('\n')]))
		message.append(red_output.getvalue())
		notify_info.emit('\n'.join(message))
	
	@classmethod
	def show_dialog(cls):
		dialog = cls()
		dialog.exec()
		
		if dialog.result() != 1:
			return
		
		cls.run_current_command()

class WidthDialog(QDialog):
	@classmethod
	def show_dialog(cls):
		dialog = cls()
		dialog.exec()

		if dialog.result():
			return
		
		value = dialog.input_width_value.value()
		expr = dialog.input_width_expression.toPlainText()

		with config:
			config['plot_width'] = value
			config['plot_widthexpression'] = expr

	def __init__(self):
		super().__init__()
		self.setWindowTitle('Width Dialog')
		self.layout = QFormLayout()
		self.setLayout(self.layout)

		self.input_width_value = QQ(QDoubleSpinBoxFullPrec, value=config['plot_width'], range=(0, None))
		self.layout.addRow('Width: ', self.input_width_value)

		expression = config['plot_widthexpression']

		self.input_width_expression = QQ(QPlainTextEdit, value=expression)
		self.layout.addRow('Expression: ', self.input_width_expression)

		buttons = QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel | QDialogButtonBox.StandardButton.Reset 
		self.button_box = QDialogButtonBox(buttons)
		self.button_box.rejected.connect(lambda: self.done(1))
		self.button_box.accepted.connect(lambda: self.done(0))
		self.button_box.button(QDialogButtonBox.StandardButton.Reset).clicked.connect(lambda _: self.reset())
		self.layout.addRow(self.button_box)

	def reset(self):
		self.input_width_value.setValue(Config.initial_values['plot_width'][0])
		self.input_width_expression.setPlainText(Config.initial_values['plot_widthexpression'][0])

	@staticmethod
	def gui_set_width(value):
		width = config['plot_width']

		if value == "+":
			width *= 3/4
		elif value == "-":
			width /= 3/4
		elif value == "++":
			width *= 1/2
		elif value == "--":
			width /= 1/2
		else:
			width *= value

		config['plot_width'] = width

class OffsetDialog(QDialog):
	@classmethod
	def show_dialog(cls):
		dialog = cls()
		dialog.exec()

		if dialog.result():
			return

		value = dialog.input_offset_value.value()
		expr = dialog.input_offset_expression.toPlainText()
		is_relative = dialog.input_offset_isrelative.isChecked()

		if not is_relative:
			ref_position = mainwindow.lwpwidget.get_current_ax().ref_position
			value -= ref_position

		with config:
			config['plot_offset'] = value
			config['plot_offsetexpression'] = expr
			config['plot_offsetisrelative'] = is_relative

	def __init__(self):
		super().__init__()
		self.setWindowTitle('Offset Dialog')
		self.layout = QFormLayout()
		self.setLayout(self.layout)

		self.input_offset_value = QQ(QDoubleSpinBoxFullPrec, value=config['plot_offset'], range=(None, None))
		self.layout.addRow('Offset: ', self.input_offset_value)

		is_relative = config['plot_offsetisrelative']
		self.input_offset_isrelative = QQ(QCheckBox, value=is_relative, text='Offset is relative')
		self.layout.addRow(self.input_offset_isrelative)

		expression = config['plot_offsetexpression']

		self.input_offset_expression = QQ(QPlainTextEdit, value=expression)
		self.layout.addRow('Expression: ', self.input_offset_expression)

		buttons = QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel | QDialogButtonBox.StandardButton.Reset 
		self.button_box = QDialogButtonBox(buttons)
		self.button_box.rejected.connect(lambda: self.done(1))
		self.button_box.accepted.connect(lambda: self.done(0))
		self.button_box.button(QDialogButtonBox.StandardButton.Reset).clicked.connect(lambda _: self.reset())
		self.layout.addRow(self.button_box)

	def reset(self):
		self.input_offset_value.setValue(Config.initial_values['plot_offset'][0])
		self.input_offset_expression.setPlainText(Config.initial_values['plot_offsetexpression'][0])
		self.input_offset_isrelative.setChecked(Config.initial_values['plot_offsetisrelative'][0])

	@staticmethod
	def gui_set_offset(value):
		offset = config['plot_offset']
		width = config['plot_width']

		if value == "+":
			offset += width/4
		elif value == "-":
			offset -= width/4
		elif value == "++":
			offset += width/2
		elif value == "--":
			offset -= width/2
		else:
			offset += value

		config['plot_offset'] = offset

class QNsDialog(QDialog):
	_instance = None
	
	def __init__(self, frequency):
		super().__init__()

		if self._instance:
			self._instance.done(0)
		self.__class__._instance = self

		noq = config["series_qns"]
		visible_qn_labels = [f'qn{ul}{n+1}' for ul in ('u', 'l') for n in range(noq)]
		self.result_ = {key: pyckett.SENTINEL for key in visible_qn_labels}


		self.setWindowTitle(f"Choose QNs for transition at {frequency}")
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

		qnslayout.setColumnStretch(noq+1, 1)
		layout.addLayout(qnslayout)

		tmp = ["dist", "x", "log y"]
		cols = tmp + visible_qn_labels
		table = QTableWidget()
		table.setColumnCount(len(cols)+1)
		table.setHorizontalHeaderLabels(["Assign"] + cols)
		table.setRowCount(0)
		table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)

		tmp_df = CatFile.get_data().copy()
		tmp_df["dist"] = tmp_df["x"] - frequency
		tmp_df["log y"] = np.log10(tmp_df["y"])
		tmp_df["absdist"] = abs(tmp_df["dist"])
		tmp_df.sort_values(by=["absdist"], inplace=True)
		tmp_df.reset_index(drop=True, inplace=True)

		for i, row in tmp_df.head(100).iterrows():
			currRowCount = table.rowCount()
			table.insertRow(currRowCount)
			for j, col in enumerate(cols):
				val = f'{row[col]:{config["flag_xformatfloat"]}}'.rstrip("0").rstrip(".")
				table.setItem(currRowCount, j+1, QTableWidgetItem(val))
			tmp_dict = {key: row[key] for key in visible_qn_labels}
			tmp_dict["xpre"] = row["x"]
			table.setCellWidget(currRowCount, 0, QQ(QPushButton, text="Assign", change=lambda x, tmp_dict=tmp_dict: self.table_save(tmp_dict)))

		table.resizeColumnsToContents()
		layout.addWidget(table)

		buttons = QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
		buttonBox = QDialogButtonBox(buttons)
		buttonBox.accepted.connect(lambda: self.selector_save(1))
		buttonBox.rejected.connect(lambda: self.done(0))
		layout.addWidget(buttonBox)

	def table_save(self, tmpd):
		self.done(1)
		self.result_.update(tmpd)

	def selector_save(self, val):
		self.done(val)
		tmp_dict = {key: sb.value() for key, sb in self.sbs.items()}
		tmp_dict["xpre"] = self.frequency
		self.result_.update(tmp_dict)

	def on_exit(self):
		self.__class__._instance = None
	
	def save(self):
		return(self.result_)

class AssignAllDialog(QDialog):
	drawplot = pyqtSignal()
	_instance = None

	@classmethod
	def show_dialog(cls, i_col=None):
		if i_col is None:
			i_col = mainwindow.lwpwidget._active_ax_index[1]
		
		dialog = cls(i_col)
		dialog.show()

	def __init__(self, i_col):
		super().__init__()
		self.setModal(False)
		self.i_col = i_col
		self.axes_to_skip = set()
		
		self.setWindowTitle("Assign All")

		self.finished.connect(self.on_exit)
		if self._instance:
			self._instance.done(0)
		self.__class__._instance = self
		

		self.init_gui()
		self.drawplot.connect(self.plot_canvas.draw_idle)
		self.update_gui()

	def update_gui_asap(self, reset_axes_to_skip=True):
		if reset_axes_to_skip:
			self.axes_to_skip = set()
		
		n_rows = config['plot_rows']
		n_qns = self.noq = config['series_qns']
		
		fit_width = config['assignall_fitwidth']
		max_fwhm = config['assignall_maxfwhm']
		max_pol_degree = config['assignall_maxdegree']
		
		peakdirection = config['fit_peakdirection']
		fitmethod = config['fit_fitmethod']
		offset = config['fit_offset']
		
		asap_axes = mainwindow.lwpwidget.lwpaxes[:, self.i_col]
		qn_labels = [f'qn{ul}{i+1}' for ul in ('u', 'l') for i in range(n_qns)]
		positions = np.zeros(asap_axes.shape)
		qns = []
		for asap_ax in asap_axes:
			qns.append(asap_ax.qns)

		# Find already assigned transitions of series
		offsets = {}
		pred_egy = {}
		for i_row, qnus in enumerate(qns):
			if qnus is None:
				continue
			
			eqy_query = ' and '.join([f'(qn{i+1} == {qn})' for i, qn in enumerate(qnus)])
			egy_offsets = ASAPAx.egy_df.query(eqy_query)['egy'].to_numpy()
			if len(egy_offsets) == 1:
				pred_egy[i_row] = egy_offsets[0]
			
			# We have to get the energy levels here as the upper levels, as this is how they are defined in the *lin format
			query = ' and '.join([f'(qnu{i+1} == {qn}) and (qnl{i+1} == 0)' for i, qn in enumerate(qnus)])
			vals = LinFile.query_c(query)['x'].to_numpy()
			
			if len(vals) == 0:
				continue
			
			if len(vals) > 1:
				qnsstring = ','.join(map(str, qnus))
				notify_warning.emit(f'Multiple assignments for level {qnsstring} found. Did not use it as reference.')
				continue

			offsets[i_row] = vals[0] - pred_egy.get(i_row, 0)
		
		if len(offsets) < 1:
			msg = 'Please assign at least a single transition of the series (two subsequent would be great).'
			notify_info.emit(msg)
			raise GUIAbortedError(msg)
		
		# Create new assignments
		for ax in self.fig.get_axes():
			self.fig.delaxes(ax)
		axs = self.fig.subplots(n_rows, sharex=True, gridspec_kw={'hspace': 0, 'wspace': 0})
		self.axs = axs = axs[::-1]
		
		self.new_assignments = {}
		for i_row, (xref, qnus, ax) in enumerate(zip(positions, qns, axs)):
			already_assigned = (i_row in offsets)
			user_flagged_to_skip = (i_row in self.axes_to_skip)
			
			if already_assigned:
				xpre = xref + offsets[i_row]
				ax.scatter(offsets[i_row], 0, color=config['color_ref'], marker="*", zorder=100)
			else:
				offset_values = np.array([[key, value] for key, value in offsets.items()])
				degree = min(len(offset_values) - 1, max_pol_degree)

				pred_pol = np.polyfit(*offset_values.T, degree)
				pred_offset = np.polyval(pred_pol, i_row)
			
				xpre = xref + pred_offset
			
			xmin, xmax = xpre - fit_width/2, xpre + fit_width/2
			
			asap_ax = asap_axes[i_row]
			corr_xs, corr_ys = asap_ax.corr_xs, asap_ax.corr_ys
			
			if corr_xs is None or corr_ys is None:
				continue
			
			minindex = corr_xs.searchsorted(xmin, side="left")
			maxindex = corr_xs.searchsorted(xmax, side="right")
			
			exp_xs, exp_ys = corr_xs[minindex:maxindex], corr_ys[minindex:maxindex]
			
			ax.plot(exp_xs, exp_ys, color=config['color_exp'])
			ax.xaxis.set_visible(False)
			ax.yaxis.set_visible(False)
			ax.margins(y=config['plot_ymargin'])
			
			if already_assigned or user_flagged_to_skip or len(exp_xs) == 0:
				continue
			
			kwargs = {'wmax': max_fwhm, 'xs_weight_factor': 0}
			fit_xs = np.linspace(xmin, xmax, 1000)
			fit_function = get_fitfunction(fitmethod, offset, kwargs=kwargs)

			try:
				xmiddle, xuncert, fit_xs, fit_ys = fit_function(exp_xs, exp_ys, peakdirection, fit_xs, )
			except Exception as E:
				notify_warning.emit(f'Error when trying to fit row {i_row}. Error reads \'{E}\'.')
				continue
			
			offsets[i_row] = xmiddle - xref
			row_qns = np.array((qnus, np.zeros_like(qnus)))
			self.new_assignments[i_row] = {'xmiddle': xmiddle, 'xuncert': xuncert, 'row_qns': row_qns, 'xref': xref, 'pred_egy': pred_egy.get(i_row, 0)}
			
			ax.plot(fit_xs - xref, fit_ys, color=config['color_fit'], lw=2, alpha=0.5)
			ax.axvline(xmiddle - xref, zorder=10, color=config['color_fit'])
		
		ax = axs[0]
		ax.xaxis.set_visible(True)
		ax.margins(x=0)
		
		self.drawplot.emit()
	
	def update_gui_llwp(self, reset_axes_to_skip=True):
		if reset_axes_to_skip:
			self.axes_to_skip = set()
		
		n_rows = config['plot_rows']
		n_qns = self.noq = config['series_qns']
		
		fit_width = config['assignall_fitwidth']
		max_fwhm = config['assignall_maxfwhm']
		max_pol_degree = config['assignall_maxdegree']
		
		peakdirection = config['fit_peakdirection']
		fitmethod = config['fit_fitmethod']
		offset = config['fit_offset']
		
		tab_widget = ReferenceSeriesWindow.instance.tab
		refwidget = tab_widget.widget(self.i_col)
		
		# @Luis: Improve for transitions case to do all available predictions 
		# Think about going giving here calc_references 100 instead of n_rows
		# Then throw out all positions that are zero
		positions, qns = refwidget.calc_references(n_rows, n_qns)
		qn_labels = [f'qn{ul}{i+1}' for ul in ('u', 'l') for i in range(n_qns)]
		
		# Find already assigned transitions of series
		offsets = {}
		for i_row, (ref_pos, row_qns) in enumerate(zip(positions, qns)):
			query = ' and '.join([f'({label} == {qn})' for qn, label in zip(row_qns.flatten(), qn_labels)])
			
			vals = LinFile.query_c(query)['x'].to_numpy()
			
			if len(vals) == 0:
				continue
			
			if len(vals) > 1:
				qnsstring = ','.join(map(str, row_qns[0])) + ' ← ' + ','.join(map(str, row_qns[1]))
				notify_warning.emit(f'Multiple assignments for transition {qnsstring} found. Did not use it as reference.')
				continue
			
			offsets[i_row] = vals[0] - ref_pos
		
		if len(offsets) < 1:
			msg = 'Please assign at least a single transition of the series (two subsequent would be great).'
			notify_info.emit(msg)
			raise GUIAbortedError(msg)
		
		
		# Create new assignments
		for ax in self.fig.get_axes():
			self.fig.delaxes(ax)
		axs = self.fig.subplots(n_rows, sharex=True, gridspec_kw={'hspace': 0, 'wspace': 0})
		self.axs = axs = axs[::-1]
		
		self.new_assignments = {}
		for i_row, (xref, row_qns, ax) in enumerate(zip(positions, qns, axs)):
			already_assigned = (i_row in offsets)
			user_flagged_to_skip = (i_row in self.axes_to_skip)
			
			if already_assigned:
				xpre = xref + offsets[i_row]
				ax.scatter(offsets[i_row], 0, color=config['color_ref'], marker="*", zorder=100)
			else:
				offset_values = np.array([[key, value] for key, value in offsets.items()])
				degree = min(len(offset_values) - 1, max_pol_degree)

				pred_pol = np.polyfit(*offset_values.T, degree)
				pred_offset = np.polyval(pred_pol, i_row)
			
				xpre = xref + pred_offset
			
			xmin, xmax = xpre - fit_width/2, xpre + fit_width/2
			df = ExpFile.get_data(xrange=(xmin, xmax)).copy()
			exp_xs, exp_ys = df['x'].to_numpy(), df['y'].to_numpy()
			
			ax.plot(exp_xs - xref, exp_ys, color=config['color_exp'])
			ax.xaxis.set_visible(False)
			ax.yaxis.set_visible(False)
			ax.margins(y=config['plot_ymargin'])
			
			if already_assigned or user_flagged_to_skip or len(exp_xs) == 0:
				continue
			
			kwargs = {'wmax': max_fwhm, 'xs_weight_factor': 0}
			fit_xs = np.linspace(xmin, xmax, 1000)
			fit_function = get_fitfunction(fitmethod, offset, kwargs=kwargs)

			try:
				xmiddle, xuncert, fit_xs, fit_ys = fit_function(exp_xs, exp_ys, peakdirection, fit_xs, )
			except Exception as E:
				notify_warning.emit(f'Error when trying to fit row {i_row}. Error reads \'{E}\'.')
				continue
			
			offsets[i_row] = xmiddle - xref
			self.new_assignments[i_row] = {'xmiddle': xmiddle, 'xuncert': xuncert, 'row_qns': row_qns, 'xref': xref}
			
			ax.plot(fit_xs - xref, fit_ys, color=config['color_fit'], lw=2, alpha=0.5)
			ax.axvline(xmiddle - xref, zorder=10, color=config['color_fit'])

		ax = axs[0]
		ax.xaxis.set_visible(True)
		ax.margins(x=0)
		
		self.drawplot.emit()

	def update_gui(self, *args, **kwargs):
		self.update_gui_llwp(*args, **kwargs)

	def init_gui(self):
		layout = QVBoxLayout(margin=True)
		self.setLayout(layout)

		self.fig = matplotlib.figure.Figure(dpi=config['plot_dpi'])
		cid = self.fig.canvas.mpl_connect('button_press_event', self.on_click)

		self.plot_canvas = FigureCanvas(self.fig)
		self.mpl_toolbar = NavigationToolbar2QT(self.plot_canvas, self)

		layout.addWidget(self.plot_canvas, 1)
		layout.addWidget(self.mpl_toolbar)

		grid_layout = QGridLayout()
		
		i_row = 0

		grid_layout.addWidget(QQ(QLabel, text='Fit Width: '), i_row, 0)
		grid_layout.addWidget(QQ(QDoubleSpinBox, 'assignall_fitwidth', range=(0, None)), i_row, 1)
		
		i_row += 1

		grid_layout.addWidget(QQ(QLabel, text='Max FWHM: '), i_row, 0)
		grid_layout.addWidget(QQ(QDoubleSpinBox, 'assignall_maxfwhm', range=(0, None)), i_row, 1)

		i_row += 1

		grid_layout.addWidget(QQ(QLabel, text='Max Pol Degree: '), i_row, 0)
		grid_layout.addWidget(QQ(QSpinBox, 'assignall_maxdegree', range=(0, None)), i_row, 1)

		layout.addLayout(grid_layout)

		buttons_layout = QHBoxLayout()

		buttons_layout.addStretch()
		buttons_layout.addWidget(QQ(QPushButton, text="Assign All", shortcut="Ctrl+Return", change=lambda x: self.assign()))
		buttons_layout.addWidget(QQ(QPushButton, text="Update", change=lambda x: self.update_gui()))
		buttons_layout.addWidget(QQ(QPushButton, text="Cancel", change=lambda x: self.close()))
		buttons_layout.addStretch()

		layout.addLayout(buttons_layout)

	def assign(self):
		self.done(0)

		xs = []
		qns = []
		errors = []
		for i_row, assignment in self.new_assignments.items():
			xs.append(assignment['xmiddle'] + assignment.get('pred_egy', 0))
			errors.append(LWPWidget._ax_class.fit_determine_uncert(assignment['xref'], assignment['xmiddle'], assignment['xuncert']))
			qns.append(assignment['row_qns'])

		qns = np.array(qns)

		new_assignments = {
			'x': xs,
			'error': errors,
			'weight': 1,
			'comment': config['fit_comment'],
			'filename': '__newassignments__',
		}

		for i in range(self.noq):
			new_assignments[f'qnu{i+1}'] = qns[:,0,i]
			new_assignments[f'qnl{i+1}'] = qns[:,1,i]
		
		for i in range(self.noq, N_QNS):
			new_assignments[f'qnu{i+1}'] = pyckett.SENTINEL
			new_assignments[f'qnl{i+1}'] = pyckett.SENTINEL

		NewAssignments.get_instance().add_rows(new_assignments)

	def on_click(self, event):
		clicked_axis = event.inaxes
		
		if not clicked_axis:
			return
		
		for i, ax in enumerate(self.axs):
			if ax == event.inaxes:
				if i in self.axes_to_skip:
					self.axes_to_skip.remove(i)
				else:
					self.axes_to_skip.add(i)
				self.update_gui(reset_axes_to_skip=False)
				return

	def on_exit(self):
		self.__class__._instance = None

##
## Additional Winodws
##
class FileWindow(EQDockWidget):
	fileaddition_requested = pyqtSignal(type, str)
	default_position = None
	available_in = ['LLWP', 'LASAP']

	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		self.setWindowTitle('Files Window')
		
		self.setAcceptDrops(True)
		
		mainwidget = QWidget()
		self.setWidget(mainwidget)
		mainlayout = QVBoxLayout()
		mainwidget.setLayout(mainlayout)
		self.tab = QTabWidget()
		mainlayout.addWidget(self.tab)
		
		self.fileclasses = {'exp': ExpFile, 'cat': CatFile, 'lin': LinFile}
		
		self.components = {key: {
			'layout': QVBoxLayout(margin=True),
			'toplayout': cls.gui_settings_general(),
			'scrollarea': QScrollArea(),
			'scrolllayout': QGridLayout(margin=True),
			'counter': AtomicCounter(),			
		} for key, cls in self.fileclasses.items()}
		
		
		for key, components in self.components.items():
			layout = components['layout']
			
			widget = QWidget()
			widget.setLayout(layout)
			scrollwidget = QWidget()
			scrollwidget.setLayout(components['scrolllayout'])
			
			scrollarea = components['scrollarea']
			layout.addLayout(components['toplayout'])
			
			scrollarea.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
			scrollarea.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
			scrollarea.setWidgetResizable(True)
			scrollarea.setWidget(scrollwidget)
			
			layout.addWidget(scrollarea)
			self.tab.addTab(widget, key.capitalize())
		
			cls = self.fileclasses[key]
			for key in cls.ids.keys():
				self.gui_add_file(cls, key)
		
		self.fileaddition_requested.connect(self.gui_add_file)

	def gui_add_file(self, cls, id):
		file = cls.ids[id]
		key = next(key for key, value in self.fileclasses.items() if value == cls)
		
		components = self.components[key]
		scrolllayout = components['scrolllayout']
		i_row = components['counter'].increase()

		row_widgets = file.gui_settings_widgets()
		scrolllayout.setRowStretch(i_row, 0)
		for i_col, widget in enumerate(row_widgets):
			scrolllayout.addWidget(widget, i_row, i_col)
		scrolllayout.setRowStretch(i_row + 1, 1)
		
	def dragEnterEvent(self, event):
		if event.mimeData().hasUrls():
			event.accept()
		else:
			event.ignore()

	def dropEvent(self, event):
		mainwindow.dropEvent(event)

class ScalingWindow(EQDockWidget):
	default_visible = False
	default_position = 2
	available_in = ['LLWP']

	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		self.setWindowTitle('Scaling Window')
		
		mainwidget = QWidget()
		self.setWidget(mainwidget)
		mainlayout = QGridLayout(margin=True)
		mainwidget.setLayout(mainlayout)

		checkbox_scale = QQ(QComboBox, "plot_yscale", options=("Per Plot", "Global", "Custom"))
		spinbox_scalemin = QQ(QDoubleSpinBox, "plot_yscale_min", range=(None, None), minWidth=120)
		spinbox_scalemax = QQ(QDoubleSpinBox, "plot_yscale_max", range=(None, None), minWidth=120)
		spinbox_scalecatfac = QQ(QDoubleSpinBox, "plot_expcat_factor", range=(0, None), minWidth=120)
		spinbox_scalecatexp = QQ(QSpinBox, "plot_expcat_exponent", range=(None, None), prefix="*10^", minWidth=120)

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
		checkbox_scale.currentTextChanged.emit(config["plot_yscale"])

		mainlayout.addWidget(QLabel("y-scale:  "), 0, 0)
		mainlayout.addWidget(checkbox_scale, 0, 1)
		mainlayout.addWidget(scaleCatLabel, 1, 0)
		mainlayout.addWidget(spinbox_scalecatfac, 1, 1)
		mainlayout.addWidget(spinbox_scalecatexp, 1, 2, 1, 2)
		mainlayout.addWidget(scaleMinLabel, 2, 0)
		mainlayout.addWidget(spinbox_scalemin, 2, 1)
		mainlayout.addWidget(scaleMaxLabel, 2, 2)
		mainlayout.addWidget(spinbox_scalemax, 2, 3)
		mainlayout.setColumnStretch(7, 10)
		mainlayout.setRowStretch(3, 1)

class ReferenceSelector(QTabWidget):
	methods = ('Transition', 'List', 'Expression')
	values_changed = pyqtSignal()
	
	def __init__(self, parent, initial_values={}):
		super().__init__(parent)
		self.parent = parent
		self.column = self.parent.tab.indexOf(self)
		self.state = {
			'method': self.methods[0],
			'check_blends': False,
			'transition': {},
			'list': {
				'qns': None,
				'xs': [],
				'i0': 0,
			},
			'expression': {
				'expression':	'',
				'N0':			0,
			},
		}
		self.state.update(initial_values)

		self.tabwidgets = {}
		for label in self.methods:
			self.tabwidgets[label] = QWidget()
			self.addTab(self.tabwidgets[label], label)

		self.setCurrentIndex(self.methods.index(self.state['method']))
		self.currentChanged.connect(lambda x: self.update_method())

		# Tab 1: Transition
		layout = QVBoxLayout(margin=True)

		text = self.state['transition'].get('file') if self.state.get('transition') else None
		text = os.path.split(text)[-1] if text else "All"
		tooltip = "Select if the first transition out of all visible files will be used or only predictions from a selected file should be used."
		
		self.activefilebutton = QQ(QPushButton, text=text, tooltip=tooltip, change=self.seriesselector_file)
		tmplayout = QHBoxLayout()
		tmplayout.addWidget(QQ(QLabel, text="File: "))
		tmplayout.addWidget(self.activefilebutton, 1)

		self.series_selector = SeriesSelector(self, self.state["transition"])
		self.series_selector.values_changed.connect(self.seriesselector_changed)
		self.seriesselector_changed(silent=True) # Call to set 'transition' field in state properly

		layout.addLayout(tmplayout)
		layout.addWidget(self.series_selector)

		# button_apply = QQ(QPushButton, text="Apply", change=lambda x: mainwindow.lwpwidget.reset_offsets())

		label_blend = QLabel("Blend Width: ")
		width_blend = QQ(QDoubleSpinBox, "series_blendwidth", minWidth=85, range=(0, None))

		button_box = QHBoxLayout()
		button_box.addWidget(label_blend)
		button_box.addWidget(width_blend)
		button_box.addStretch(1)

		layout.addLayout(button_box)
		layout.addStretch(1)
		self.tabwidgets["Transition"].setLayout(layout)


		# Tab 2: List
		layout = QVBoxLayout(margin=True)

		self.listtable = QQ(QTableWidget, rowCount=0, columnCount=2, move=(0, 0))
		self.listtable.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
		self.listtable.setHorizontalHeaderLabels(["#", "Frequency"])
		layout.addWidget(self.listtable)
		
		if len(self.state['list']['xs']):
			self.list_load(values=self.state["list"]["xs"])

		button_open = QQ(QToolButton, text="Open List", change=lambda x: self.list_load(dialog=False))
		button_write = QQ(QToolButton, text="Write List", change=lambda x: self.list_load(dialog=True))

		label_startat = QLabel("Start at Index: ")
		spinbox_startat = QQ(QSpinBox, value=self.state["list"]["i0"], range=(0, None), singlestep=1, change=self.list_indexupdate)
		self.spinbox_startat = spinbox_startat
		# button_apply = QQ(QPushButton, text="Apply", change=lambda x: mainwindow.lwpwidget.reset_offsets())

		hbox = QHBoxLayout()
		[ hbox.addWidget(label_startat), hbox.addWidget(spinbox_startat), hbox.addStretch(1)]
		layout.addLayout(hbox)

		hbox = QHBoxLayout()
		[hbox.addWidget(button_open), hbox.addWidget(button_write), hbox.addStretch(1)]
		layout.addLayout(hbox)

		self.tabwidgets["List"].setLayout(layout)

		# Tab 3: Expression
		layout = QVBoxLayout(margin=True)

		placeholder = "Enter expression here\ne.g. for a linear molecule with B=4000 use: (N+N0)*4000*2"
		self.input_expr = QQ(QPlainTextEdit, value=self.state["expression"]["expression"], placeholder=placeholder, change=lambda: self.state["expression"].__setitem__("expression", self.input_expr.toPlainText()))		
		# button_apply = QQ(QPushButton, text="Apply", change=lambda x: mainwindow.lwpwidget.reset_offsets())
		button_update = QQ(QPushButton, text="Update", change=lambda x: self.expression_update(type='expresion'))

		label_N0 = QQ(QLabel, text="N0: ")
		self.input_N0 = QQ(QSpinBox, value=self.state["expression"]["N0"], range=(0, None), change=lambda x: self.expression_update(type='index'))

		layout.addWidget(self.input_expr)
		hbox = QHBoxLayout()
		[hbox.addWidget(label_N0),
		hbox.addWidget(self.input_N0), hbox.addStretch(1), hbox.addWidget(button_update)]
		layout.addLayout(hbox)

		self.tabwidgets["Expression"].setLayout(layout)

	def update_method(self):
		current_tab = self.currentIndex()
		self.state['method'] = self.methods[current_tab]
		self.changed()

	def seriesselector_file(self):
		options = ["All", *CatFile.ids.keys()]
		current_file = self.state['transition'].get("file")
		if current_file in options:
			current = options.index(current_file)
		else:
			current = 0
		item, ok = QInputDialog.getItem(self, "Choose File", f"Limit this series to the following file:", options, current=current, editable=False)
		if not (ok and item):
			return

		if item == "All":
			item, itemshort = None, "All"
		else:
			itemshort = os.path.split(item)[-1]

		self.state['transition']['file'] = item
		self.series_selector.state['file'] = item
		
		self.activefilebutton.setText(itemshort)
		self.changed()

	def seriesselector_changed(self, silent=False):
		self.state['transition'] = self.series_selector.state.copy()
		if not silent:
			self.changed()

	def list_indexupdate(self):
		self.state['list']['i0'] = self.spinbox_startat.value()
		self.changed()

	def list_load(self, values=None, dialog=False):
		if values is not None:
			xs = values
			self.state["list"]["xs"] = xs

		elif dialog:
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
					notify_warning.emit(f"Could not convert the string '{x}' to a numerical value.")

			self.state["list"].update({
				"qns":		None,
				"xs":		xs,
			})
			self.spinbox_startat.setValue(0)

		else:
			fname, _ = QFileDialog.getOpenFileName(self, 'Open Positions List',)
			if not fname:
				return
			xs = []
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
							notify_warning.emit(f"Could not convert the string '{x}' to a numerical value.")

			self.state["list"].update({
				"qns":		None,
				"xs":		xs,
			})
			self.spinbox_startat.setValue(0)


		table = self.listtable
		table.setRowCount(0)
		i=0
		for i, x in enumerate(xs):
			currRowCount = table.rowCount()
			table.insertRow(currRowCount)
			table.setItem(currRowCount, 0, QTableWidgetItem(f"{i}"))
			table.setItem(currRowCount, 1, QTableWidgetItem(f"{x:{config['flag_xformatfloat']}}"))
		self.changed()

	def expression_update(self, type=None):
		if type == 'index':
			self.state['expression']['N0'] = self.input_N0.value()
		else:
			expression = self.input_expr.toPlainText()
			expression = expression.replace('\r', '').replace('\n', '')
			self.state['expression']['expresssion'] = expression

		self.changed()

	def seriesselector_getqns(self, n_qns):
		state = self.state['transition']
		tmp = {
			'qnus': state['qnus'],
			'qnls': state['qnls'],
			'diff': state['diff'] if state['use_diff'] else state['incr'],
		}

		tmp = {key: np.array(value[:n_qns]) for key, value in tmp.items()}
		return(tmp)

	def calc_references(self, n_rows, n_qns):
		method = self.state['method']
		shape = (n_rows, 2, n_qns)
		qns = np.zeros(shape, dtype=np.int64)
		positions = np.zeros(n_rows, dtype=np.float64)
		self.state['check_blends'] = False

		if method == 'Transition':
			state = self.state['transition']
			tmp = self.seriesselector_getqns(n_qns)
			qnus, qnls, diffs = tmp['qnus'], tmp['qnls'], tmp['diff']

			# Prefilter df to all transitins belonging to series
			conditions, conditions_incr = [], []
			normalizing_value = None
			for i, qnu, qnl, diff in zip(range(n_qns), qnus, qnls, diffs):
				diff = int(diff)
				if diff:
					if normalizing_value is None:
						normalizing_value = qnu // diff
					conditions_incr.append(f"((qnu{i+1} - {qnu-normalizing_value*diff})/{diff})")
					conditions_incr.append(f"((qnl{i+1} - {qnl-normalizing_value*diff})/{diff})")
				else:
					conditions.append(f"(qnu{i+1} == {qnu})")
					conditions.append(f"(qnl{i+1} == {qnl})")

			if len(conditions_incr):
				conditions.append(" == ".join(conditions_incr))
			
			# Add condition for limiting to specific file
			file_to_limit_data_to = state.get('file')
			if file_to_limit_data_to:
				filepath = file_to_limit_data_to.replace('\\', '\\\\')
				conditions.append(f'(filename == "{filepath}")')
			
			conditions = " and ".join(conditions)  if conditions else '(visible)'
			cat_df = CatFile.query_c(conditions)

			# Select the desired transitions
			for i_row in range(n_rows):
				tmp_qnus, tmp_qnls = qnus + i_row * diffs, qnls + i_row * diffs
				qns[i_row] = (tmp_qnus, tmp_qnls)

				cond_upper = [f"(qnu{i+1} == {qn})" for i, qn in enumerate(tmp_qnus)]
				cond_lower = [f"(qnl{i+1} == {qn})" for i, qn in enumerate(tmp_qnls)]
				condition  = " & ".join(cond_upper + cond_lower)

				vals = cat_df.query(condition)["x"].to_numpy()
				val = vals[0] if len(vals) else 0
				positions[i_row] = val
			
			self.state['check_blends'] = True
		elif method == "List":
			state = self.state['list']
			i0 = state['i0']
			xs = state['xs']
			tmp_qns = state['qns']

			if xs is not None and i0 < len(xs):
				imax = min(len(positions), len(xs)-i0)
			else:
				imax = 0
			
			positions[:imax] = xs[i0:imax+i0]

			if tmp_qns is not None and len(tmp_qns) and i0 < len(tmp_qns):
				qns[:imax] = qns[i0:imax+i0]
		
		else: # Expression
			state = self.state['expression']
			N0 = state['N0']
			expression = state['expression']


			if ';'in expression: # expression gives QNs
				expressions = [expr.strip() for expr in expression.split(';') if expr.strip()]
				if len(expressions) % 2 != 0:
					raise CustomError('Error with the provided expression. Number of subexpressions has to be even.')
				n_qns = len(expressions) // 2
				shape = (n_rows, 2, n_qns)
				qns = np.zeros(shape, dtype=np.int64)

				try:
					expressions = [compile(expr, '', 'eval') for expr in expressions]
					rows_qns = [ [eval(expr, {"N": i, "N0": N0}) for expr in expressions] for i in range(n_rows)]
				except Exception as E:
					notify_error.emit('The provided expression could not be evaluated.')
					raise

				# Prefilter predictions to used series
				diff_qns = np.diff(np.array(rows_qns), axis=0)
				diffs = [np.unique(col) for col in diff_qns.T]
				qnus, qnls = np.split(np.array(rows_qns[0]), 2)

				conditions, conditions_incr = [], []
				normalizing_value = None
				for i, qnu, qnl, diff in zip(range(n_qns), qnus, qnls, diffs):
					if len(diff) != 1:
						continue
					diff = int(diff[0])
					if diff:
						if normalizing_value is None:
							normalizing_value = qnu // diff
						conditions_incr.append(f"((qnu{i+1} - {qnu-normalizing_value*diff})/{diff})")
						conditions_incr.append(f"((qnl{i+1} - {qnl-normalizing_value*diff})/{diff})")
					else:
						conditions.append(f"(qnu{i+1} == {qnu})")
						conditions.append(f"(qnl{i+1} == {qnl})")

				if len(conditions_incr):
					conditions.append(" == ".join(conditions_incr))

				conditions = " and ".join(conditions) if conditions else '(visible)'
				cat_df = CatFile.query_c(conditions)

				# Get the exact predictions from prefiltered dataframe
				qn_labels = [f'qn{ul}{i+1}' for ul in ('u', 'l') for i in range(n_qns)]
				for i_row, row_qns in enumerate(rows_qns):
					query = " and ".join([f'({label} == {qn})' for qn, label in zip(row_qns, qn_labels) ])
					vals = cat_df.query(query)["x"].to_numpy()
					val = vals[0] if len(vals) else 0
					positions[i_row] = val
					qns[i_row] = np.split(np.array(row_qns), 2)

				self.state['check_blends'] = False
			else: # expression gives xs
				if expression:
					try:
						expression = compile(expression, '', 'eval')
						positions = [eval(expression, {"N": i, "N0": N0}) for i in range(n_rows)]
					except Exception as E:
						notify_error.emit('The provided expression could not be evaluated.')
		
		return(positions, qns)

	def changed(self):
		self.values_changed.emit()
		mainwindow.lwpwidget.set_data()

	def contextMenuEvent(self, event):
		menu = QMenu(self)
		get_positions_action = menu.addAction('Copy Reference Positions')
		get_qns_action = menu.addAction('Copy Reference QNs')
		fit_all_action = menu.addAction('Fit all')

		action = menu.exec(self.mapToGlobal(event.pos()))
		if action == get_positions_action:
			n_rows = config['plot_rows']
			n_qns = config['series_qns']
			positions, _ = self.calc_references(n_rows, n_qns)
			QApplication.clipboard().setText('\n'.join(map(str, positions[::-1])))
		elif action == get_qns_action:
			n_rows = config['plot_rows']
			n_qns = config['series_qns']
			_, qns = self.calc_references(n_rows, n_qns)

			output_string = []
			for qnus, qnls in qns[::-1]:
				tmp_string = ','.join([f'{qn:3.0f}' for qn in qnus]) + ' ← ' + ','.join([f'{qn:3.0f}' for qn in qnls])
				output_string.append(tmp_string)

			output_string = '\n'.join(output_string)
			QApplication.clipboard().setText(output_string)
		elif action == fit_all_action:
			i_col = self.parent.tab.indexOf(self)
			AssignAllDialog.show_dialog(i_col)

class SeriesSelector(QWidget):
	values_changed = pyqtSignal()

	def __init__(self, parent, initial_values={}):
		super().__init__(parent)
		self.n_qns = N_QNS
		self.parent = parent
		self.updating = False
		self.state = {
			'qnus': (1, 0, 1, 0, 0, 0),
			'qnls': (0, 0, 0, 0, 0, 0),
			'incr': (True, False, True, False, False, False),
			'diff': (1, 0, 1, 0, 0, 0),
			'use_diff': False,
		}
		self.state.update(initial_values)

		layout = QGridLayout()

		create_qn = lambda: QQ(QSpinBox, minWidth=60, maxWidth=60, range=(None, None), visible=False,
							singlestep=1, change=lambda x: self.changed())

		self.qnus = [create_qn() for _ in range(self.n_qns)]
		self.qnls = [create_qn() for _ in range(self.n_qns)]

		create_widget = lambda widget, kwargs: QQ(widget, minWidth=40, maxWidth=40, visible=False,
							change=lambda x: self.changed(), **kwargs)

		self.incr = [create_widget(QCheckBox, {'text':'Inc'}) for _ in range(self.n_qns)]
		self.diff = [create_widget(QSpinBox, {'range': (None, None), 'singlestep': 1, }) for _ in range(self.n_qns)]

		self.incqns = QQ(QPushButton, text="Inc", change=lambda x: self.incdecqns(+1), width=40)
		self.decqns = QQ(QPushButton, text="Dec", change=lambda x: self.incdecqns(-1), width=40)

		self.togglediff = QQ(QToolButton, text="⇆", change=lambda x: self.change_incr_mode(), width=40)
		
		for i, widget in enumerate(self.qnus + self.qnls):
			layout.addWidget(widget, i//self.n_qns, i%self.n_qns)

		for i, incr, diff in zip(range(N_QNS), self.incr, self.diff):
			tmp = QHBoxLayout()
			tmp.addWidget(incr)
			tmp.addWidget(diff)

			layout.addLayout(tmp, 4, i)
			layout.setColumnStretch(i, 100)

		for i in range(self.n_qns):
			layout.setColumnStretch(i, 100)
		
		layout.addWidget(self.togglediff, 4, self.n_qns, 1, 1)

		layout.addWidget(self.incqns, 0, self.n_qns, 1, 2)
		layout.addWidget(self.decqns, 1, self.n_qns, 1, 2)

		layout.setRowStretch(6, 10)
		layout.setColumnStretch(self.n_qns+2, 1)

		self.layout = layout
		self.setLayout(layout)
		self.set_state()

		config.register_widget("series_qns", self.togglediff, self.set_state)
		config.register_widget("flag_showseriesarrows", self.togglediff, self.change_arrows)
		self.change_arrows()

	def set_state(self):
		self.updating = True
		state = self.state

		n_qns = config["series_qns"]
		are_visible = [True] * n_qns +[False] * (self.n_qns - n_qns)

		for qnu, value, is_visible in zip(self.qnus, state["qnus"], are_visible):
			qnu.setValue(value)
			qnu.setVisible(is_visible)
		for qnl, value, is_visible in zip(self.qnls, state["qnls"], are_visible):
			qnl.setValue(value)
			qnl.setVisible(is_visible)
		for widget, value, is_visible in zip(self.incr, state["incr"], are_visible):
			widget.setChecked(value)
			widget.setVisible(is_visible and not state['use_diff'])
		for widget, value, is_visible in zip(self.diff, state["diff"], are_visible):
			widget.setValue(value)
			widget.setVisible(is_visible and state['use_diff'])
		self.updating = False
		self.changed()

	def change_incr_mode(self):
		self.state['use_diff'] = not self.state['use_diff']
		self.set_state()

	def incdecqns(self, dir):
		self.updating = True
		incr_values = (x.value() for x in self.diff) if self.state['use_diff'] else (x.isChecked() for x in self.incr)

		for qnu, qnl, incr in zip(self.qnus, self.qnls, incr_values):
			qnu.setValue(qnu.value()+dir*incr)
			qnl.setValue(qnl.value()+dir*incr)
		self.updating = False
		self.changed()

	def changed(self):
		if self.updating:
			return
	
		self.state["qnus"] = [x.value() for x in self.qnus]
		self.state["qnls"] = [x.value() for x in self.qnls]
		self.state["incr"] = [x.isChecked() for x in self.incr]
		self.state["diff"] = [x.value() for x in self.diff]

		self.values_changed.emit()

	def wheelEvent(self, event):
		steps = event.angleDelta().y() // 120
		self.incdecqns(steps)
	
	def change_arrows(self):
		show_arrows = config['flag_showseriesarrows']
		width = 60 if show_arrows else 40
		button_symbols = QAbstractSpinBox.ButtonSymbols.UpDownArrows if show_arrows else QAbstractSpinBox.ButtonSymbols.NoButtons

		for widget in self.qnus + self.qnls + self.diff:
			widget.setMinimumWidth(width)
			widget.setMaximumWidth(width)
			widget.setButtonSymbols(button_symbols)


class ReferenceSeriesWindow(EQDockWidget):
	default_visible = True
	available_in = ['LLWP']

	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		self.setWindowTitle("Reference Series")

		widget = QGroupBox()
		layout = QVBoxLayout()
		self.setWidget(widget)
		widget.setLayout(layout)

		self.tab = QTabWidget()
		self.tab_order = None
		layout.addWidget(self.tab)

		self.tab.setTabsClosable(True)
		self.tab.setMovable(True)
		self.tab.setDocumentMode(True)

		self.tab.setTabBarAutoHide(True)
		self.tab.setCornerWidget(QQ(QToolButton, text="Dupl.", tooltip="Duplicate current tab", change=self.duplicate_tab), Qt.Corner.TopRightCorner)
		self.tab.tabCloseRequested.connect(self.close_tab)
		self.tab.tabBarDoubleClicked.connect(self.renameoradd_tab)
		self.tab.setCurrentIndex(config["series_currenttab"])

		self.tab_order = [self.tab.widget(i) for i in range(self.tab.count())]
		self.set_state(config["series_references"])
		if not self.tab.count():
			self.add_tab()

		mainwindow.lwpwidget.plotscreated.connect(self.update_number_of_references)
		self.tab.currentChanged.connect(self.check_order)

	def close_tab(self, index, check_order=True):
		if self.tab.count() == 1:
			return
		tab = self.tab.widget(index)
		tab.deleteLater()
		self.tab.removeTab(index)
		if check_order:
			self.check_order()

	def renameoradd_tab(self, index):
		if index == -1:
			self.add_tab()
		elif self.tab.widget(index) != 0:
			text, ok = QInputDialog().getText(self, "Tab Name","Enter the Tabs Name:")
			if ok and text:
				self.tab.setTabText(index, text)
				self.tab.widget(index).state["title"] = text

	def add_tab(self, init_values={}, check_order=True):
		title = init_values.get("title", "Series")
		tmp = ReferenceSelector(self, init_values)
		tmp.values_changed.connect(self.changed)
		self.tab.addTab(tmp, title)
		if check_order:
			self.check_order()

	def duplicate_tab(self, *args, **kwargs):
		values = self.tab.widget(self.tab.currentIndex()).state.copy()
		self.add_tab(values, *args, **kwargs)

	def update_number_of_references(self):
		numberoftabs = self.tab.count()
		diff = config["plot_cols"] - numberoftabs
		locked = config["flag_syncreferencestocolumns"]

		if diff > 0:
			for i in range(diff):
				self.duplicate_tab(check_order=False)
		elif diff < 0 and locked:
			for i in range(-diff):
				self.close_tab(numberoftabs - i - 1, check_order=False)
		
	def set_state(self, values):
		for i in range(self.tab.count()):
			tab = self.tab.widget(i)
			tab.deleteLater()
		self.tab.clear()
		
		for tabdata in values:
			self.add_tab(tabdata)

	def get_state(self):
		values = []
		for i in range(self.tab.count()):
			tmp = self.tab.widget(i).state
			values.append(tmp)
		return(values)
	
	@drawplot_decorator.d
	def check_order(self, x=None):
		if self.tab_order is None:
			return

		if x is None:
			x = self.tab.currentIndex()
		config['series_currenttab'] = x
		new_tab_order = [self.tab.widget(i) for i in range(self.tab.count())]
		
		if new_tab_order == self.tab_order:
			return

		self.changed()
		mainwindow.lwpwidget.set_data()

	
	def redo_references(self):
		self.tab_order = []
		self.check_order()

	def changed(self):
		tmp = self.get_state()
		config["series_references"] = tmp

class LogWindow(EQDockWidget):
	available_in = ['LLWP', 'LASAP']

	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		self.setWindowTitle("Log")

		mainwidget = QGroupBox()
		layout = QVBoxLayout()
		self.setWidget(mainwidget)
		mainwidget.setLayout(layout)

		self.log_area = QTextEdit()
		self.log_area.setReadOnly(True)
		self.log_area.setMinimumHeight(50)

		notify_error.connect(lambda x: self.writelog(x, style='error'))
		notify_warning.connect(lambda x: self.writelog(x, style='warning'))
		notify_info.connect(lambda x: self.writelog(x, style='info'))

		layout.addWidget(self.log_area)

	def writelog(self, text, style=None):
		separator = "<br/>"
		tmp = self.log_area.toHtml()
		tmp = tmp.split(separator)
		if len(tmp)-1 > config["flag_logmaxrows"]:
			self.log_area.setHtml(separator.join(tmp[-config["flag_logmaxrows"]:]))

		time_str = time.strftime("%H:%M", time.localtime())
		

		style_string = {
			'error': "<span style='color:#ff0000;'>ERROR</span>: ",
			'warning': "<span style='color:#eda711;'>WARNING</span>: ",
			'info': "<span style='color:#0096FF;'>INFO</span>: ",
		}.get(style, '')

		text = f"{time_str}: {style_string}{text}{separator}"
		self.log_area.append(text)
		sb = self.log_area.verticalScrollBar()
		sb.setValue(sb.maximum())

class CatTableModel(QAbstractTableModel):
	def __init__(self, headers, get_df, table):
		super().__init__()
		self.columns = list(headers.keys())
		self.headers = list(headers.values())
		self.get_df = get_df
		self.get_df_visible = lambda: get_df().loc[:, self.columns]
		self.table = table

		full_df = self.get_df()
		view_df = self.get_df_visible()
		col_pos_full = {key: i for i, key in enumerate(full_df.columns)}
		col_pos_view = {key: i for i, key in enumerate(view_df.columns)}
		self.col_trans_dict = {i_col: col_pos_full[key] for key, i_col in col_pos_view.items()}

		self.resizing_programmatically = False
		self.columns_skip_autoresize = set()
		self.table.horizontalHeader().sectionResized.connect(self.exclude_from_autoresize)

		config.register("series_qns", self.update_columns_visibility)

	def data(self, index, role):
		if role not in (Qt.ItemDataRole.DisplayRole, Qt.ItemDataRole.EditRole):
			return

		df = self.get_df_visible()
		value = df.iloc[index.row(), index.column()]

		if isinstance(value, str):
			return(value)
		elif isinstance(value, (np.integer, int)):
			if value == pyckett.SENTINEL:
				return("")
			else:
				if role == Qt.ItemDataRole.EditRole:
					return(str(value))
				else:
					return(f"{{:{config['flag_tableformatint']}}}".format(value))
		elif np.isnan(value):
			return("")
		else: #float
			if role == Qt.ItemDataRole.EditRole:
				return(str(value))
			else:
				return(f"{{:{config['flag_tableformatfloat']}}}".format(value))

	def rowCount(self, index):
		return(self.get_df_visible().shape[0])

	def columnCount(self, index):
		return(self.get_df_visible().shape[1])

	def headerData(self, section, orientation, role):
		if not role == Qt.ItemDataRole.DisplayRole:
			return
		
		if orientation == Qt.Orientation.Horizontal:
			return str(self.headers[section])

		if orientation == Qt.Orientation.Vertical:
			df = self.get_df_visible()
			if section >= len(df.index):
				return ""
			return str(df.index[section])

	def flags(self, index):
		return(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable | Qt.ItemFlag.ItemIsEditable)

	def update(self):
		self.layoutChanged.emit()

	def setData(self, index, value, role):
		if not index.isValid() or role != Qt.ItemDataRole.EditRole:
			return False

		df = self.get_df_visible()

		i_row = index.row()
		if not 0 <= i_row < df.shape[0]:
			return False
		i_col = index.column()
		if not 0 <= i_col < df.shape[1]:
			return False

		dtype = df[df.columns[i_col]].dtypes
		
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

		orig_df = self.get_df()
		orig_i_col = self.col_trans_dict[i_col]
		orig_df.iloc[i_row, orig_i_col] = value
		self.dataChanged.emit(index, index)
		self.resize_columns()
		return True

	def exclude_from_autoresize(self, index):
		if not  self.resizing_programmatically:
			self.columns_skip_autoresize.add(index)
	
	def resize_columns(self, reset=False):
		self.resizing_programmatically = True
		if reset:
			self.columns_skip_autoresize = set()
		for i_col in range(len(self.columns)):
			if i_col not in self.columns_skip_autoresize:
				self.table.resizeColumnToContents(i_col)
		self.resizing_programmatically = False

	def update_columns_visibility(self):
		qns = config['series_qns']
		for i in range(N_QNS):
			self.table.setColumnHidden(i,   i>=qns)
			self.table.setColumnHidden(i+N_QNS, i>=qns)

class NewAssignmentsWindow(EQDockWidget):
	default_visible = True
	default_position = 2
	available_in = ['LLWP', 'LASAP']

	def __init__(self):
		super().__init__()
		self.setWindowTitle("New Assignments")

		mainwidget = QGroupBox()
		layout = QVBoxLayout()
		self.setWidget(mainwidget)
		mainwidget.setLayout(layout)

		tooltip_append = "Append to file if checked or overwrite content if unchecked"
		tooltip_uncert = "Positive Values are absolute Values, -1 for obs-calc, -2 for dialog, -3 for StdDev from Fit"

		new_assignments = self.new_assignments = NewAssignments.get_instance()

		widgets = self.widgets = {
			'save': QQ(QToolButton, text='Save', change=lambda x: new_assignments.save_gui()),
			'delete': QQ(QToolButton, text='×', change=lambda x: self.delete()),
			'delete_all': QQ(QToolButton, text='Del All', change=lambda x: self.delete_all()),
			'add_row': QQ(QToolButton, text='+', change=lambda x: addemptyrow_inplace(new_assignments.get_new_assignments_df(), self.model)),
			'resize': QQ(QToolButton, text='Resize', change=lambda x: self.model.resize_columns(reset=True)),
			# 'append': QQ(QCheckBox, 'flag_appendonsave', text='Append', tooltip=tooltip_append),
			'error_label': QQ(QLabel, text='Default Uncertainty: ', tooltip=tooltip_uncert),
			'error': QQ(QDoubleSpinBox, 'fit_uncertainty', range=(-3, None), minWidth=120, singlestep=config["fit_uncertaintystep"], tooltip=tooltip_uncert),
		}

		tmp = { f'qn{UL.lower()}{i+1}': f'{UL}{i+1}' for UL in 'UL' for i in range(N_QNS) }
		headers = { **tmp, 'x': 'Freq', 'error': 'Unc', 'weight': 'Weight', 'comment': 'Comment'}

		self.table = QTableView()
		self.model = CatTableModel(headers, new_assignments.get_new_assignments_df, self.table)
		self.table.setModel(self.model)
		self.model.update_columns_visibility()
		self.model.resize_columns()

		buttonsBox = QHBoxLayout()
		for key in ('delete', 'delete_all', 'add_row', None, 'save', 'resize'):
			buttonsBox.addWidget(widgets[key]) if key is not None else buttonsBox.addStretch(1)

		layout.addLayout(buttonsBox)
		layout.addWidget(self.table)
		buttonsBox = QHBoxLayout()
		
		for key in ('error_label', None, 'error'):
			buttonsBox.addWidget(widgets[key]) if key is not None else buttonsBox.addStretch(2)
		layout.addLayout(buttonsBox)

	def scroll_bottom(self):
		self.table.selectRow(len(self.model.get_df())-1)
		self.table.scrollToBottom()

	def delete_all(self):
		self.delete(delete_all=True)

	def delete(self, delete_all=False):
		df = self.new_assignments.get_new_assignments_df()
		df.reset_index(drop=True, inplace=True)

		if delete_all:
			df.drop(df.index, inplace=True)
		else:
			selected = [index.row() for index in self.table.selectionModel().selectedRows()]
			df.drop(selected, inplace=True)

		df.reset_index(inplace=True, drop=True)
		
		self.model.update()
		self.new_assignments.load_file()
		mainwindow.lwpwidget.set_data()

class ConvolutionWindow(EQDockWidget):
	default_visible = False
	default_position = None
	available_in = ['LLWP',]


	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		self.setWindowTitle('Lineshape Window')

		mainwidget = QWidget()
		self.setWidget(mainwidget)

		self.general_layout = QGridLayout(margin=True)
		mainwidget.setLayout(self.general_layout)
		row_i = 0

		options = ('Off', 'Gauss', 'Lorentz', 'Voigt', 'Custom')
		self.selection_box = QQ(QComboBox, 'convolution_function', options=options)
		self.selection_box_label = QQ(QLabel, text='Function: ')
		self.general_layout.addWidget(self.selection_box_label, row_i, 0)
		self.general_layout.addWidget(self.selection_box, row_i, 1)

		row_i += 1

		self.stepwidth_input = QQ(QDoubleSpinBox, 'convolution_stepwidth', range=(0, None))
		self.stepwidth_label = QQ(QLabel, text='Stepwidth: ')
		self.general_layout.addWidget(self.stepwidth_label, row_i, 0)
		self.general_layout.addWidget(self.stepwidth_input, row_i, 1)

		row_i += 1

		self.kernelwidth_input = QQ(QDoubleSpinBox, 'convolution_kernelwidth', range=(0, None))
		self.kernelwidth_label = QQ(QLabel, text='Kernel Width: ')
		self.general_layout.addWidget(self.kernelwidth_label, row_i, 0)
		self.general_layout.addWidget(self.kernelwidth_input, row_i, 1)

		row_i += 1

		self.derivative_input = QQ(QSpinBox, 'convolution_derivative', range=(0, 2))
		self.derivative_label = QQ(QLabel, text='Derivative: ')
		self.general_layout.addWidget(self.derivative_label, row_i, 0)
		self.general_layout.addWidget(self.derivative_input, row_i, 1)

		row_i += 1

		self.amplitude_input = QQ(QDoubleSpinBox, 'convolution_amplitude', range=(None, None))
		self.amplitude_label = QQ(QLabel, text='Amplitude: ')
		self.general_layout.addWidget(self.amplitude_label, row_i, 0)
		self.general_layout.addWidget(self.amplitude_input, row_i, 1)

		row_i += 1

		self.width_gauss_input = QQ(QDoubleSpinBox, 'convolution_widthgauss', range=(None, None))
		self.width_gauss_label = QQ(QLabel, text='Width Gauss: ')
		self.general_layout.addWidget(self.width_gauss_label, row_i, 0)
		self.general_layout.addWidget(self.width_gauss_input, row_i, 1)

		row_i += 1

		self.width_lorentz_input = QQ(QDoubleSpinBox, 'convolution_widthlorentz', range=(None, None))
		self.width_lorentz_label = QQ(QLabel, text='Width Lorentz: ')
		self.general_layout.addWidget(self.width_lorentz_label, row_i, 0)
		self.general_layout.addWidget(self.width_lorentz_input, row_i, 1)

		row_i += 1

		self.custom_kernel_input = QQ(QPlainTextEdit)
		self.custom_kernel_label = QQ(QLabel, text='Custom Kernel: ')
		self.general_layout.addWidget(self.custom_kernel_label, row_i, 0)
		self.general_layout.addWidget(self.custom_kernel_input, row_i, 1)

		row_i += 1

		self.update_button = QQ(QPushButton, text='Update', change=lambda _: self.on_change())
		self.general_layout.addWidget(self.update_button, row_i, 0, 1, 2)

		row_i += 1

		self.general_layout.setRowStretch(row_i, 1)

		keys = ('convolution_function', 'convolution_stepwidth', 'convolution_kernelwidth', 
		  'convolution_derivative', 'convolution_amplitude', 'convolution_widthgauss', 
		  'convolution_widthlorentz')
		config.register(keys, self.on_change)
		config.register('convolution_function', self.update_gui)

		self.update_gui()

	def update_gui(self):
		function = config['convolution_function'] 

		self.stepwidth_input.setVisible(function != 'Off')
		self.stepwidth_label.setVisible(function != 'Off')

		self.kernelwidth_input.setVisible(function not in ('Off', 'Custom'))
		self.kernelwidth_label.setVisible(function not in ('Off', 'Custom'))

		self.derivative_input.setVisible(function != 'Off')
		self.derivative_label.setVisible(function != 'Off')

		self.amplitude_input.setVisible(function != 'Off')
		self.amplitude_label.setVisible(function != 'Off')

		self.width_gauss_input.setVisible(function in ('Gauss', 'Voigt'))
		self.width_gauss_label.setVisible(function in ('Gauss', 'Voigt'))

		self.width_lorentz_input.setVisible(function in ('Lorentz', 'Voigt'))
		self.width_lorentz_label.setVisible(function in ('Lorentz', 'Voigt'))

		self.custom_kernel_input.setVisible(function == 'Custom')
		self.custom_kernel_label.setVisible(function == 'Custom')

		self.update_button.setVisible(function == 'Custom')

	def on_change(self):
		function = config['convolution_function']

		if function == 'Off':
			kernel = []
		elif function in ('Gauss', 'Lorentz', 'Voigt'):
			stepwidth = config['convolution_stepwidth']
			kernelwidth = config['convolution_kernelwidth']
			derivative = config['convolution_derivative']
			amplitude = config['convolution_amplitude']

			number_of_points = kernelwidth // (2 * stepwidth)
			kernel_xs = np.arange(-number_of_points, +number_of_points+1) * stepwidth

			if function == 'Gauss':
				width = config['convolution_widthgauss']
				kernel = lineshape('Gauss', derivative, kernel_xs, 0, amplitude, width)
			elif function == 'Lorentz':
				width = config['convolution_widthlorentz']
				kernel = lineshape('Lorentz', derivative, kernel_xs, 0, amplitude, width)
			else: # Voigt
				widths = config['convolution_widthgauss'], config['convolution_widthlorentz']
				kernel = lineshape('Voigt', derivative, kernel_xs, 0, amplitude, *widths)
		else: # Custom Kernel
			tmp = self.custom_kernel_input.toPlainText().strip()
			if not tmp:
				kernel = []
			else:
				tmp = tmp.split(',')
				kernel = [float(x) for x in tmp]

		config['convolution_kernel'] = list(kernel)
		mainwindow.lwpwidget.set_data()

class CloseByLinesWindow(EQDockWidget):
	cursor_changed = pyqtSignal(float, float)
	default_visible = False
	default_position = None
	available_in = ['LLWP',]


	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		self.setWindowTitle('Close-By Lines')

		mainwidget = QWidget()
		self.setWidget(mainwidget)

		self.layout = QVBoxLayout(margin=True)
		mainwidget.setLayout(self.layout)

		self.text_is_frozen_label = QQ(QLabel)

		toplayout = QHBoxLayout()
		toplayout.addWidget(QQ(QLabel, text='Cutoff Width: '))
		toplayout.addWidget(QQ(QDoubleSpinBox, 'plot_hovercutoff', range=(0, None)))
		toplayout.addStretch(1)
		toplayout.addWidget(self.text_is_frozen_label)

		self.layout.addLayout(toplayout)

		self.text_area = QTextEdit()
		self.text_area.setReadOnly(True)
		self.text_area.setMinimumHeight(50)

		self.layout.addWidget(self.text_area)
		self.cursor_changed.connect(self.on_cursor_change)
		
		self._text_is_frozen = False

		tmp = QShortcut('Ctrl+Shift+F', mainwindow)
		tmp.activated.connect(self.toggle_text_is_frozen)
		tmp.setContext(Qt.ShortcutContext.ApplicationShortcut)
	
	def toggle_text_is_frozen(self):
		self.text_is_frozen = not self.text_is_frozen
	
	def on_cursor_change(self, x, y):
		if not self.isVisible():
			return
		
		if self.text_is_frozen:
			return
		self.update_closeby_lines(x)

	def update_closeby_lines(self, xcenter):
		cutoff = config['plot_hovercutoff']
		xrange = xcenter - cutoff, xcenter + cutoff

		cat_df = CatFile.get_data(xrange=xrange)
		lin_df = LinFile.get_data(xrange=xrange)

		dataframes = {'cat': cat_df, 'lin': lin_df}
		textrows = []
		noq = config['series_qns']

		for type, df in dataframes.items():
			if not len(df):
				continue

			for row in df.to_dict(orient="records"):
				qnus = [row[f'qnu{i+1}'] for i in range(noq)]
				qnls = [row[f'qnl{i+1}'] for i in range(noq)]

				qnus_string = ''.join([f'{qn:3.0f}' for qn in qnus if qn != pyckett.SENTINEL])
				qnls_string = ''.join([f'{qn:3.0f}' for qn in qnls if qn != pyckett.SENTINEL])

				qnstring = f'{qnus_string} ← {qnls_string}'

				if type == 'cat':
					ylog = np.log10(row['y'])
					text = config['closebylines_catfstring'].format(qns=qnstring, ylog=ylog, **row)
				else: # lin
					text = config['closebylines_linfstring'].format(qns=qnstring, **row)

				textrows.append((row['x'], text))

		fmt = config['flag_tableformatfloat']
		textrows.append((xcenter, f'<span id="anchor">CURSOR POSITION: {xcenter:{fmt}}</span>'))

		textrows = sorted(textrows, key=lambda x: x[0])
		total_text = '<br>'.join([x[1] for x in textrows])
		total_text = f'<pre>{total_text}</pre>'
		self.text_area.setHtml(total_text)
		self.text_area.scrollToAnchor(f'anchor')

	@property
	def text_is_frozen(self):
		return(self._text_is_frozen)
	
	@text_is_frozen.setter
	def text_is_frozen(self, value):
		self._text_is_frozen = value
		text = 'Text is frozen' if self._text_is_frozen else 'Text is unfrozen'
		self.text_is_frozen_label.setText(text)

class ConfigWindow(EQDockWidget):
	default_visible = False
	default_position = None
	available_in = ['LLWP', 'LASAP']

	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		self.setWindowTitle('Config')

		vbox = QVBoxLayout(margin=True)
		scrollarea = QScrollArea()
		widget = QWidget()
		layout = QGridLayout(margin=True)

		self.updating = True

		scrollarea.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
		scrollarea.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
		scrollarea.setWidgetResizable(True)

		tmp_layout = QHBoxLayout()
		tmp_layout.addWidget(QQ(QPushButton, text='Save as default', change=lambda: config.save()))
		completer = QCompleter(config.keys())
		completer.setCaseSensitivity(Qt.CaseSensitivity.CaseInsensitive)
		tmp_layout.addWidget(QQ(QLineEdit, placeholder="Search", completer=completer, change=lambda x: self.search(x)))
		tmp_layout.addStretch(1)

		vbox.addLayout(tmp_layout)
		self.widgets = {}

		i = 1
		for key, value in config.items():
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

		mainwidget = QWidget()
		self.setWidget(mainwidget)
		mainwidget.setLayout(vbox)

		self.updating = False
		self.visibilityChanged.connect(self.on_visibility_change)

	def on_visibility_change(self, is_visible):
		if is_visible:
			self.timer = QTimer(self)
			self.timer.timeout.connect(self.get_values)
			self.timer.start(200)
		else:
			self.timer.stop()

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
			value = config[key]
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
		converter = Config.initial_values.get(key)
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
			config[key] = value
			oklab.setText("Good")
		except Exception as E:
			oklab.setText("Bad")


class CreditsWindow(EQDockWidget):
	default_visible = False
	default_position = None
	available_in = ['LLWP', 'LASAP']

	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		self.setWindowTitle("Credits")

		global CREDITSSTRING
		self.setWidget(QQ(QLabel, text=CREDITSSTRING, align=Qt.AlignmentFlag.AlignCenter, wordwrap=True, minHeight=400, minWidth=500))


class ResidualsWindow(EQDockWidget):
	default_visible = False
	default_position = None
	available_in = ['LLWP',]

	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		self.setWindowTitle('Residuals')

		self.fig = matplotlib.figure.Figure(dpi=config['plot_dpi'])
		self.plot_canvas = FigureCanvas(self.fig)

		self.ax = self.fig.subplots()

		self.points = self.ax.scatter([], [], color=[], marker=".")
		self.annot = self.ax.annotate("", xy=(0,0), xytext=(5,5), textcoords="offset points", color="black", ha="center", va="bottom", bbox=dict(boxstyle="round", fc="w"))
		self.annot.set_visible(False)
		self.fig.canvas.mpl_connect("motion_notify_event", self.on_hover)
		self.fig.canvas.mpl_connect("button_press_event", lambda x: self.on_hover(x, True))

		self.mpl_toolbar = NavigationToolbar2QT(self.plot_canvas, self)
		self.df = None

		layout = QVBoxLayout(margin=True)
		mainwidget = QWidget()
		self.setWidget(mainwidget)
		mainwidget.setLayout(layout)
		layout.addWidget(self.plot_canvas, 6)
		hlayout = QHBoxLayout()
		hlayout.addWidget(QLabel("x-axis: "))
		hlayout.addWidget(QQ(QLineEdit, "residuals_xvariable", placeholder="Choose the x-variable, e.g. x_lin"))
		hlayout.addWidget(QLabel("y-axis: "))
		hlayout.addWidget(QQ(QLineEdit, "residuals_yvariable", placeholder="Choose the y-variable, e.g. x_lin-x_cat"))
		hlayout.addWidget(QQ(QCheckBox, "residuals_blends", text="Blends"))
		hlayout.addWidget(QQ(QCheckBox, "residuals_autoscale", text="Autoscale on Update"))
		layout.addLayout(hlayout)
		layout.addWidget(self.mpl_toolbar)
		layout.addWidget(QQ(QPlainTextEdit, "residuals_query", maxHeight=60, placeholder="Query text to filter shown lines. Use qnu1, ..., qnu6 and qnl1, ..., qnl6 for the quantum numbers. Other possible values are x_lin, x_cat, error_lin, error_cat, filename_lin, filename_cat, y, degfreed, elower, usd, tag, qnfmt, weight, and comment."))
		layout.addWidget(QQ(QPlainTextEdit, "residuals_colorinput", maxHeight=60, placeholder="Enter custom color and query to color specific lines differently. E.g. enter '#ff0000; qnu1 < 20' to color all transitions with the first upper quantum number below 20 red."))

		buttonslayout = QHBoxLayout()
		layout.addLayout(buttonslayout)
		buttonslayout.addStretch(1)
		self.update_button = QQ(QPushButton, text="Update", change=self.plot_residuals)
		buttonslayout.addWidget(self.update_button)
		buttonslayout.addWidget(QQ(QPushButton, text="Save", change=self.save_residuals))
		buttonslayout.addStretch(1)

	def get_residuals(self):
		lin_df = LinFile.get_data()
		cat_df = CatFile.get_data()

		noq = config["series_qns"]
		self.noq = noq
		qns_visible = [f"qn{ul}{n+1}" for ul in ("u", "l") for n in range(noq)]
		df = pd.merge(lin_df, cat_df, how="inner", on=qns_visible)
		df.rename(columns={"x_x": "x_lin", "x_y": "x_cat", "error_x": "error_lin", "error_y": "error_cat", "filename_x": "filename_lin", "filename_y": "filename_cat"}, inplace=True)

		query = config["residuals_query"].strip()
		if query:
			df.query(query, inplace=True)
		df.reset_index(drop=True, inplace=True)

		if config["residuals_blends"]:
			mask = df["weight"] != 0
			df_view = df[mask]
			tmp_dict = df_view.groupby(df_view.x_lin).apply(lambda x: np.average(x.x_cat, weights=x.weight)).to_dict()
			df["x_cat"] = df["x_lin"].map(lambda x: tmp_dict.get(x, x))

		df["obs_calc"] = df["x_lin"] - df["x_cat"]
		return(df)

	def plot_residuals(self):
		self.update_button.setDisabled(True)

		try:
			thread = self.plot_residuals_core()
			thread.wait()
		except Exception as E:
			notify_warning.emit("There was an error in your Residuals window input")
		finally:
			self.fig.canvas.draw_idle()
			self.update_button.setDisabled(False)

	@QThread.threaded_d
	def plot_residuals_core(self, thread=None):
		message = []

		df = self.get_residuals()
		df["color"] = config["residuals_defaultcolor"]
		message.append(f"Found {len(df)} entries matching your query.")

		colorquerytext = config["residuals_colorinput"].split("\n")
		for row in colorquerytext:
			if row.strip():
				color, query = row.split(";")
				indices = df.query(query).index
				df.loc[indices, "color"] = color
				message.append(f"{len(indices)} lines are colored in <span style='color:{color};'>{color}</span>.")

		self.df = df

		yvariable = config["residuals_yvariable"].strip() or "obs_calc"
		xvariable = config["residuals_xvariable"].strip() or "x_lin"
		ys = df.eval(yvariable).to_numpy()
		xs = df.eval(xvariable).to_numpy()
		colors = df["color"].to_numpy()
		tuples = list(zip(xs,ys))
		tuples = tuples if len(tuples)!=0 else [[None,None]]
		self.points.set_offsets(tuples)
		self.points.set_color(colors)
		if len(xs) and config["residuals_autoscale"]:
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
			self.ax.set_ylim(y_range[0]-config["plot_ymargin"]*(y_range[1]-y_range[0]), y_range[1]+config["plot_ymargin"]*(y_range[1]-y_range[0]))
		notify_info.emit("<br/>".join(message))

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

			# clicked and transition(s) under cursor
			if click and noq:
				action_to_perform = 'showinplot'
				
				is_right_click = (event.button ==3)
				if is_right_click:
					menu = QMenu()
					showinplot_action = menu.addAction('Show in plot')
					deletefromfile_action = menu.addAction('Delete from file')
					
					action = menu.exec(QCursor.pos())
					if action == showinplot_action:
						action_to_perform = 'showinplot'
					elif action == deletefromfile_action:
						action_to_perform = 'deletefromfile'

				if action_to_perform == 'showinplot':
					tab_widget = ReferenceSeriesWindow.instance.tab
					refwidget = tab_widget.widget(config['series_currenttab'])
					refwidget.setCurrentIndex(0)
					seriesselector = refwidget.series_selector
					
					current_state = seriesselector.state
					current_state['qnus'][:noq] = qnus
					current_state['qnls'][:noq] = qnls
					config['series_qns'] = noq

					seriesselector.set_state()

					mainwindow.lwpwidget.set_data()
					mainwindow.show()
					mainwindow.raise_()
					mainwindow.activateWindow()
				elif action_to_perform == 'deletefromfile':
					unique_lin_fnames = tmp_transitions['filename_lin'].unique()
					for lin_fname in unique_lin_fnames:
						if lin_fname == "__newassignments__":
							new_assignments = NewAssignments.get_instance()
							lin = new_assignments.get_new_assignments_df().copy()
						else:
							lin = pyckett.lin_to_df(lin_fname, sort=False)
						
						transitions_to_delete = tmp_transitions.query('filename_lin == @lin_fname')
						columns = {x: x.replace('_lin', '').replace('_x', '') for x in transitions_to_delete.columns if (x.endswith('_lin') or x.endswith('_x'))}
						transitions_to_delete = transitions_to_delete.rename(columns=columns)
						transitions_to_delete = transitions_to_delete[list(pyckett.lin_dtypes.keys())]

						lin = pd.merge(lin, transitions_to_delete, how='left', indicator=True)
						lin = lin[lin['_merge'] == 'left_only']
						lin = lin.drop('_merge', axis='columns').reset_index(drop=True)

						if lin_fname == '__newassignments__':
							new_assignments.new_assignments_df = lin
							new_assignments.load_file()

							new_assignments_window = NewAssignmentsWindow.instance
							new_assignments_window.model.update()
							new_assignments_window.model.resize_columns()
						else:
							with open(lin_fname, 'w+') as file:
								file.write(pyckett.df_to_lin(lin))


class BlendedLinesWindow(EQDockWidget):
	default_visible = False
	default_position = None
	available_in = ['LLWP',]

	set_indicator_text = pyqtSignal(str)
	fill_table_requested = pyqtSignal()

	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		self.setWindowTitle('Blended Lines')

		self.peaks = []
		self.fit_values = None
		self.cid = None

		mainwidget = QWidget()
		self.setWidget(mainwidget)
		layout = QVBoxLayout(margin=True)
		mainwidget.setLayout(layout)
		
		class CustomPlotWidget(PlotWidget):
			update_plot_requested = pyqtSignal()
			
			def gui(self):
				super().gui()
				self.fit_line = self.ax.plot([], [], color = config['blendedlines_color_total'])[0]
				self.plot_parts = []
				self.plot_parts_lock = threading.Lock()
				self.update_plot_requested.connect(super().update_plot)
					
			def update_plot(self):
				self.parent.fit_peaks(self)
		
		self.plot_widget = CustomPlotWidget(parent=self)

		layout.addWidget(self.plot_widget, 2)

		tmplayout = QGridLayout()
		row_id = 0
		tmplayout.addWidget(QQ(QLabel, text="Function: "), row_id, 0)
		tmplayout.addWidget(QQ(QComboBox, "blendedlines_lineshape", items=("Gauss", "Lorentz", "Voigt"), minWidth=120), row_id, 1)

		tmplayout.addWidget(QQ(QLabel, text="Transparency: "), row_id, 2)
		tmplayout.addWidget(QQ(QDoubleSpinBox, "blendedlines_transparency", range=(0, 1), minWidth=120, singlestep=0.1), row_id, 3)

		row_id += 1
		tmplayout.addWidget(QQ(QLabel, text="Derivative: "), row_id, 0)
		tmplayout.addWidget(QQ(QSpinBox, "blendedlines_derivative", range=(0, 2), minWidth=120), row_id, 1)

		row_id += 1
		tmplayout.addWidget(QQ(QLabel, text="Max FWHM: "), row_id, 0)
		tmplayout.addWidget(QQ(QDoubleSpinBox, "blendedlines_maxfwhm", range=(0, None), minWidth=120), row_id, 1)

		tmplayout.addWidget(QQ(QCheckBox, "blendedlines_fixedwidth", text="All Same Width"), row_id, 2, 1, 2)

		row_id += 1
		tmplayout.addWidget(QQ(QLabel, text="Baseline Rank: "), row_id, 0)
		tmplayout.addWidget(QQ(QSpinBox, "blendedlines_polynom", range=(-1, None), minWidth=120), row_id, 1)

		tmplayout.addWidget(QQ(QCheckBox, "blendedlines_showbaseline", text="Show Baseline"), row_id, 2, 1, 2)
		tmplayout.setColumnStretch(4, 1)

		layout.addLayout(tmplayout)

		self.label = QQ(QLabel, text="Ready", textFormat=Qt.TextFormat.RichText)
		self.set_indicator_text.connect(self.label.setText)
		tmp_layout = QHBoxLayout()
		layout.addLayout(tmp_layout)
		row = (
		  QQ(QPushButton, text="Del All", change=lambda x: self.del_peak(-1)),
		  QQ(QPushButton, text="Update", change=lambda x: self.plot_widget.update_plot()),
		  QQ(QPushButton, text="Save", change=lambda x: self.save_values()),
		  self.label,
		)

		for widget in row:
			tmp_layout.addWidget(widget)
		tmp_layout.addStretch()

		self.table = QTableWidget()
		self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
		self.table.setMinimumHeight(50)
		layout.addWidget(self.table, 1)

		self.cid = self.plot_widget.plot_canvas.mpl_connect("button_press_event", lambda event: self.add_peak(event))
		self.fill_table_requested.connect(self.fill_table)
		self.plot_widget.from_current_plot()

	@QThread.threaded_d
	def fit_peaks(self, plot_widget, thread=None):
		self.set_indicator_text.emit("<span style='color:#eda711;'>Working ...</span>")
		
		peaks = self.peaks.copy()
		profile = config["blendedlines_lineshape"]
		derivative = config["blendedlines_derivative"]
		polynomrank = config["blendedlines_polynom"]+1
		fixedwidth = config["blendedlines_fixedwidth"]
		amplitude_direction = config["fit_peakdirection"]
		now = 2 if profile == "Voigt" else 1
		noa = 2 + now * (not fixedwidth)

		fitfunction = lambda x, *args, fixedwidth=fixedwidth: self.fitfunction(x, profile, derivative, polynomrank, *args, fixedwidth=fixedwidth)
		fitfunction_no_baseline = lambda x, *args, fixedwidth=fixedwidth: self.fitfunction(x, profile, derivative, 0, *args, fixedwidth=fixedwidth)

		with plot_widget.plot_parts_lock:
			for part in plot_widget.plot_parts:
				part.remove()
			plot_widget.plot_parts = []

		thread.earlyreturn()

		xrange = xmin, xmax = plot_widget.xrange
		xwidth = xmax - xmin
		xcenter = (xmin + xmax) / 2

		df_exp = ExpFile.get_data(xrange=xrange).copy()
		exp_xs = df_exp["x"].to_numpy()
		exp_ys = df_exp["y"].to_numpy()

		thread.earlyreturn()

		xs = []
		ys = []
		ws = []
		exp_mean = 0
		if len(exp_ys) and (polynomrank + len(peaks)):
			if config["blendedlines_autopositionpeaks"]:
				exp_mean = exp_ys.mean()

			yptp = 4*(np.amax(exp_ys)-np.amin(exp_ys))
			w0 = xwidth
			wmax = config["blendedlines_maxfwhm"] or w0
			amp_min, amp_max = -3*yptp, 3*yptp
			if amplitude_direction < 0:
				amp_max = 0
			if amplitude_direction > 0:
				amp_min = 0
			
			for peak in peaks:
				x, y, x_rel = peak

				if not xmin < x < xmax:
					x = xcenter + x_rel
					if not xmin < x < xmax:
						x = xcenter
					y = yptp * np.sign(amplitude_direction)
				elif not amp_min < y < amp_max:
					y = yptp * np.sign(amplitude_direction)

				xs.append((x, *xrange))
				ys.append((y, amp_min, amp_max))
				ws.append((0, 0, wmax))

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
			res_xs = np.linspace(xmin, xmax, config["blendedlines_xpoints"])
			res_ys = fitfunction(res_xs, *popt)
			res_exp_ys = fitfunction(exp_xs, *popt)
		else:
			popt = [0]*(noa*len(peaks)+polynomrank+2*fixedwidth)
			perr = [0]*(noa*len(peaks)+polynomrank+2*fixedwidth)
			res_xs = np.linspace(xmin, xmax, config["blendedlines_xpoints"])
			res_ys = res_xs*0
			res_exp_ys = exp_xs*0

		thread.earlyreturn()


		ax = plot_widget.ax

		plot_widget.fit_line.set_data(res_xs, res_ys)
		plot_widget.fit_line.set_color(config["blendedlines_color_total"])

		opt_param = []
		err_param = []

		for i in range(len(peaks)):
			tmp_params = list(popt[i*noa: (i+1)*noa])
			tmp_errors = list(perr[i*noa: (i+1)*noa])

			if fixedwidth:
				tmp_params.extend(popt[-(polynomrank+now):len(popt)-polynomrank])
				tmp_errors.extend(perr[-(polynomrank+now):len(popt)-polynomrank])

			tmp_ys = fitfunction_no_baseline(res_xs, *tmp_params)
			tmp_ys += exp_mean
			
			with plot_widget.plot_parts_lock:
				plot_widget.plot_parts.append(ax.plot(res_xs, tmp_ys, color=config["blendedlines_color_total"], alpha=config["blendedlines_transparency"])[0])

			opt_param.append( tmp_params )
			err_param.append( tmp_errors )

		with plot_widget.plot_parts_lock:
			plot_widget.plot_parts.append(ax.scatter([x[0] for x in opt_param], [x[1] + exp_mean for x in opt_param], color=config["blendedlines_color_points"]))

		if polynomrank > 0 and config["blendedlines_showbaseline"]:
			baseline_args = popt[-polynomrank:]
			with plot_widget.plot_parts_lock:
				plot_widget.plot_parts.append(ax.plot(res_xs, np.polyval(baseline_args, res_xs-xcenter), color=config["blendedlines_color_baseline"])[0])
		else:
			baseline_args = []

		rms_ys = np.sqrt(np.sum((exp_ys - res_exp_ys)**2) / len(exp_ys)) if len(exp_ys) else 0
		self.params = opt_param, err_param, profile, derivative, noa, now, xcenter, baseline_args, rms_ys

		thread.earlyreturn()

		self.fill_table_requested.emit()
		self.set_indicator_text.emit("Ready")
		plot_widget.update_plot_requested.emit()
	
	def add_peak(self, event):
		x, y = event.xdata, event.ydata
		is_left_mouse_button = (event.button == 1)
		if not (x and y and is_left_mouse_button):
			return
		x = np.float64(x)
		y = np.float64(y)
		xmin, xmax = self.plot_widget.xrange
		xcenter = (xmax + xmin) / 2
		x_rel = x - xcenter

		self.peaks.append((x, y, x_rel))
		self.peaks.sort(key=lambda x: x[0])
		self.plot_widget.update_plot()

	def del_peak(self, i=None):
		if i == -1:
			self.peaks = []
		elif i is not None and isinstance(i, (int, float, np.integer, np.float64)):
			if i in range(len(self.peaks)):
				del self.peaks[int(i)]
		else:
			if len(self.peaks) != 0:
				self.peaks.pop()
		self.plot_widget.update_plot()

	def fitfunction(self, x, fun, der, bpr, *ps, fixedwidth=False):
		noa = 2 + (1+(fun == "Voigt")) * (not fixedwidth)
		now = 2 if fun == "Voigt" else 1
		res_ys = []
		param_peaks = ps

		# Baseline Polynom
		if bpr > 0:
			param_baseline = ps[-bpr:]
			param_peaks    = ps[:-bpr]
			xcenter = sum(self.plot_widget.xrange)/2

			res_ys.append(np.polyval(param_baseline, x-xcenter))

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
			table.setCellWidget(currRowCount, 0, QQ(QPushButton, text="Assign", change=lambda x, xpos=x, error=x_error: self.assign(xpos, error)))
			table.setItem(currRowCount, 1, QTableWidgetItem(f'{x:{config["flag_xformatfloat"]}}'))
			table.setItem(currRowCount, 2, QTableWidgetItem(f'{y:{config["flag_xformatfloat"]}}'))
			table.setItem(currRowCount, 3, QTableWidgetItem(f'{wg:{config["flag_xformatfloat"]}}'))
			table.setItem(currRowCount, 4, QTableWidgetItem(f'{wl:{config["flag_xformatfloat"]}}'))
			table.setCellWidget(currRowCount, 5, QQ(QPushButton, text="Assign other QNs", change=lambda x, xpos=x, error=x_error: self.assign_other_qns(xpos, error)))
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
			notify_info.emit(f"Saved the fit to the file {filename}.")
		else:
			notify_info.emit(f"No fit values to be saved.")

	
	def assign(self, x, error):
		index = self.plot_widget.index
		try:
			lwpax = mainwindow.lwpwidget.lwpaxes[index]
		except IndexError as E:
			message = 'The reference ax is not existing anymore. Assigning is therefore not possible.'
			notify_error.emit(message)
			raise GUIAbortedError(message)
		
		new_assignment = {'x': x, 'error': lwpax.fit_determine_uncert(lwpax.ref_position, x, error), 'xpre': lwpax.ref_position}
		new_assignment.update(lwpax.create_qns_dict(complete=True))
		new_assignment.update({'weight': 1, 'comment': config['fit_comment'], 'filename': '__newassignments__'})
		
		if lwpax.check_blends(new_assignment):
			return()
		NewAssignments.get_instance().add_row(new_assignment)
	
	def assign_other_qns(self, x, error):
		index = self.plot_widget.index
		try:
			lwpax = mainwindow.lwpwidget.lwpaxes[index]
		except IndexError as E:
			message = 'The reference ax is not existing anymore. Assigning is therefore not possible.'
			notify_error.emit(message)
			raise GUIAbortedError(message)
		
		new_assignment = {'x': x, 'error': lwpax.fit_determine_uncert(lwpax.ref_position, x, error), 'xpre': lwpax.ref_position}
		new_assignment.update(lwpax.create_qns_dict(complete=True))
		new_assignment.update({'weight': 1, 'comment': config['fit_comment'], 'filename': '__newassignments__'})
		
		dialog = QNsDialog(x)
		dialog.exec()

		if dialog.result() != 1:
			return
		
		new_assignment.update(dialog.save())
		if lwpax.check_blends(new_assignment):
			return()
		NewAssignments.get_instance().add_row(new_assignment)
	

	def activateWindow(self):
		if hasattr(self, 'plot_widget'):
			self.plot_widget.from_current_plot()
		super().activateWindow()

	def closeEvent(self, *args, **kwargs):
		return super().closeEvent(*args, **kwargs)

class ReportWindow(EQDockWidget):
	default_visible = False
	default_position = None
	available_in = ['LLWP',]

	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		self.setWindowTitle("Report Window")

		layout = QVBoxLayout(margin=True)
		mainwidget = QWidget()
		mainwidget.setLayout(layout)
		self.setWidget(mainwidget)

		layout.addWidget(QQ(QPlainTextEdit, "report_query", maxHeight=40, placeholder="Query text to filter assignments. Use qnu1, ..., qnu6 and qnl1, ..., qnl6 for the quantum numbers. Other possible values are x, error, weight, comment, and filename."))
		layout.addWidget(QQ(QCheckBox, "report_blends", text="Blends"))
		layout.addWidget(QQ(QPushButton, text="Create Report", change=self.create_report))
		self.reportfield = QQ(QTextEdit, readonly=True)
		self.reportfield.setFontFamily("Courier")
		layout.addWidget(self.reportfield)

	def create_report(self):
		results = {}
		report = []

		lin_df = LinFile.get_data()
		cat_df = CatFile.get_data()
		noq = config["series_qns"]
		qns_visible = [f"qn{ul}{n+1}" for ul in ("u", "l") for n in range(noq)]
		blends = config["report_blends"]

		query = config["report_query"].strip()
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

class SeriesfinderWindow(EQDockWidget):
	default_visible = False
	default_position = None
	available_in = ['LLWP',]

	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		self.setWindowTitle("Series Finder")

		layout = QVBoxLayout(margin=True)
		mainwidget = QWidget()
		mainwidget.setLayout(layout)
		self.setWidget(mainwidget)
		
		self.messageLabel = QQ(QLabel, wordwrap=True, hidden=True)

		self.outputTable = QTableWidget()
		self.outputTable.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)

		vertLayout = QHBoxLayout()
		leftLayout = QGridLayout()
		rightLayout = QVBoxLayout()
		layout.addWidget(QQ(QLabel, wordwrap=True, text="Series Finder allows to find the strongest (unassigned) predicted transitions."))

		rightLayout.addWidget(QQ(QLabel, text="Allowed Transitions"))
		rightLayout.addWidget(QQ(QCheckBox, "seriesfinder_atype", text="a-type"))
		rightLayout.addWidget(QQ(QCheckBox, "seriesfinder_btype", text="b-type"))
		rightLayout.addWidget(QQ(QCheckBox, "seriesfinder_ctype", text="c-type"))
		rightLayout.addWidget(QQ(QCheckBox, "seriesfinder_xtype", text="x-type"))
		rightLayout.addStretch(1)

		leftLayout.addWidget(QQ(QLabel, text="Start Frequency: "), 1, 0)
		leftLayout.addWidget(QQ(QLineEdit, "seriesfinder_start"), 1, 1)
		leftLayout.addWidget(QQ(QLabel, text="Stop Frequency: "), 2, 0, 1, 1)
		leftLayout.addWidget(QQ(QLineEdit, "seriesfinder_stop"), 2, 1, 1, 1)
		leftLayout.addWidget(QQ(QLabel, text="Number of Results: "), 3, 0, 1, 1)
		leftLayout.addWidget(QQ(QSpinBox, "seriesfinder_results", range=(0, None)), 3, 1, 1, 1)
		leftLayout.addWidget(QQ(QLabel, text="Additional Condition: "), 4, 0, 1, 1)
		leftLayout.addWidget(QQ(QLineEdit, "seriesfinder_condition", placeholder="Hover for tooltip", tooltip="Use qnu1, ..., qnu6, qnl1, ..., qnl6 for the quantum numbers. Additionally, x, error, y, degfreed, elower, usd, tag, qnfmt, and filename are allowed."), 4, 1, 1, 1)
		leftLayout.addWidget(QQ(QCheckBox, "seriesfinder_onlyunassigned", text = "Only unassigned Lines"), 5, 1, 1, 1)
		leftLayout.addWidget(QQ(QPushButton, text="Run", change=lambda x: self.run()), 6, 1, 1, 1)

		leftLayout.setColumnStretch(1, 1)

		vertLayout.addLayout(leftLayout)
		vertLayout.addLayout(rightLayout)
		layout.addLayout(vertLayout)
		layout.addWidget(self.messageLabel)
		layout.addWidget(self.outputTable, 1)

	def run(self):
		nor = config["seriesfinder_results"]

		condition = []
		tmp_min = config["seriesfinder_start"]
		tmp_max = config["seriesfinder_stop"]
		addCondition = config["seriesfinder_condition"]

		if addCondition.strip():
			condition.append(addCondition)
		if tmp_min:
			condition.append(f"{tmp_min} <= x")
		if tmp_max:
			condition.append(f"x <= {tmp_max}")

		tmp_condition = []
		if config["seriesfinder_atype"]:
			tmp_condition.append(f"(abs(qnu2-qnl2) % 2 == 0 and abs(qnu3-qnl3) % 2 == 1)")
		if config["seriesfinder_btype"]:
			tmp_condition.append(f"(abs(qnu2-qnl2) % 2 == 1 and abs(qnu3-qnl3) % 2 == 1)")
		if config["seriesfinder_ctype"]:
			tmp_condition.append(f"(abs(qnu2-qnl2) % 2 == 1 and abs(qnu3-qnl3) % 2 == 0)")
		if config["seriesfinder_xtype"]:
			tmp_condition.append(f"(abs(qnu2-qnl2) % 2 == 0 and abs(qnu3-qnl3) % 2 == 0)")
		if tmp_condition:
			condition.append(" or ".join(tmp_condition))
		condition = " and ".join([f"({x})" for x  in condition])

		tmp_cat_df = CatFile.get_data().copy()

		if condition:
			try:
				tmp_cat_df.query(condition, inplace=True)
			except Exception as E:
				notify_warning.emit(f"There is a syntax error in your condition: {str(E)}")
				return

		self.noq = noq = config["series_qns"]

		qns_visible = [f"qn{ul}{n+1}" for ul in ("u", "l") for n in range(noq)]
		qns_invisible = [f"qn{ul}{n+1}" for ul in ("u", "l") for n in range(noq, N_QNS)]

		if config["seriesfinder_onlyunassigned"]:
			tmp_lin_df = LinFile.get_data().copy()
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
			table.setItem(currRowCount, 1, QTableWidgetItem(f'{row["y"]:{config["flag_xformatfloat"]}}'))
			table.setItem(currRowCount, 2, QTableWidgetItem(f'{row["x"]:{config["flag_xformatfloat"]}}'))

			for i, column in enumerate(qns_visible):
				tmp = row[column]
				if tmp == pyckett.SENTINEL:
					tmp = ""
				else:
					tmp = f'{tmp:g}'
				table.setItem(currRowCount, i+3, QTableWidgetItem(tmp))

		self.outputTable.resizeColumnsToContents()

	def startHere(self, row):
		qnus = [int(row[f"qnu{i+1}"]) for i in range(self.noq)]
		qnls = [int(row[f"qnl{i+1}"]) for i in range(self.noq)]
		
		tab_widget = ReferenceSeriesWindow.instance.tab
		refwidget = tab_widget.widget(config['series_currenttab'])
		refwidget.setCurrentIndex(0)
		seriesselector = refwidget.series_selector
		
		current_state = seriesselector.state
		current_state['qnus'][:self.noq] = qnus
		current_state['qnls'][:self.noq] = qnls
		config['series_qns'] = self.noq

		seriesselector.set_state()

		mainwindow.lwpwidget.set_data()
		mainwindow.show()
		mainwindow.raise_()
		mainwindow.activateWindow()

class EnergyLevelsWindow(EQDockWidget):
	default_visible = False
	default_position = None
	available_in = ['LLWP',]

	plotting_started = pyqtSignal()
	plotting_finished = pyqtSignal()

	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		self.setWindowTitle("Energy Levels")
		self.setAcceptDrops(True)

		layout = QVBoxLayout(margin=True)
		mainwidget = QWidget()
		mainwidget.setLayout(layout)
		self.setWidget(mainwidget)

		self.fig = matplotlib.figure.Figure(dpi=config["plot_dpi"])
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

		layout.addWidget(self.plot_canvas, 6)
		layout.addWidget(self.mpl_toolbar)
		hlayout = QHBoxLayout()
		hlayout.addWidget(QQ(QPushButton, text="Open", change=lambda x: self.load_file()))
		self.file_label = QQ(QLabel, text="No File loaded")
		hlayout.addWidget(self.file_label)
		hlayout.addWidget(QLabel("x-axis: "))
		hlayout.addWidget(QQ(QLineEdit, "energylevels_xvariable", placeholder="Choose the x-variable, e.g. qn1"))
		hlayout.addWidget(QLabel("y-axis: "))
		hlayout.addWidget(QQ(QLineEdit, "energylevels_yvariable", placeholder="Choose the x-variable, e.g. egy"))
		hlayout.addWidget(QQ(QCheckBox, "energylevels_autoscale", text="Autoscale on Update"))
		layout.addLayout(hlayout)
		layout.addWidget(QQ(QPlainTextEdit, "energylevels_query", maxHeight=40, placeholder="Query text to filter shown levels. Use qn1, ..., qn6 for the quantum numbers. Other possible values are iblk, indx, egy, err, pmix, and we."))
		layout.addWidget(QQ(QPlainTextEdit, "energylevels_colorinput", maxHeight=40, placeholder="Enter custom color and query to color specific lines differently. E.g. enter '#ff0000; qn1 < 20' to color all levels with the first quantum number below 20 red."))

		buttonslayout = QHBoxLayout()
		layout.addLayout(buttonslayout)
		buttonslayout.addStretch(1)
		self.update_button = QQ(QPushButton, text="Update", change=lambda: self.plot_energylevels())
		buttonslayout.addWidget(self.update_button)
		buttonslayout.addStretch(1)
		
		self.plotting_started.connect(lambda: self.update_button.setDisabled(True))
		self.plotting_finished.connect(lambda: self.update_button.setDisabled(False))
		self.plotting_finished.connect(lambda: self.fig.canvas.draw_idle())

	def load_file(self, fname=None):
		if fname is None:
			fname = QFileDialog.getOpenFileName(None, 'Choose Egy File to load',"")[0]
		if fname:
			self.dataframe = pyckett.egy_to_df(fname)
			self.dataframe_filtered = self.dataframe
			self.fname = fname
			self.file_label.setText(os.path.split(fname)[1])
			self.plot_energylevels()

	@QThread.threaded_d
	def plot_energylevels(self, thread=None):
		self.plotting_started.emit()

		self.noq = config["series_qns"]
		try:
			if self.dataframe is None:
				return
			df = self.dataframe
			query = config["energylevels_query"]
			if query:
				df = df.query(query).copy()

			df.loc[:, "color"] = config["energylevels_defaultcolor"]
			colorquerytext = config["energylevels_colorinput"].split("\n")
			for row in colorquerytext:
				if row.strip():
					color, query = row.split(";")
					df.loc[df.query(query).index, "color"] = color

			self.dataframe_filtered = df

			xvariable = config["energylevels_xvariable"].strip() or "qn1"
			yvariable = config["energylevels_yvariable"].strip() or "egy"
			xs = df.eval(xvariable).to_numpy()
			ys = df.eval(yvariable).to_numpy()
			colors = df["color"].to_numpy()
			tuples = list(zip(xs,ys))
			tuples = tuples if len(tuples)!=0 else [[None,None]]
			self.points.set_offsets(tuples)
			self.points.set_color(colors)
			if len(xs) and config["energylevels_autoscale"]:
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
				self.ax.set_ylim(y_range[0]-config["plot_ymargin"]*(y_range[1]-y_range[0]), y_range[1]+config["plot_ymargin"]*(y_range[1]-y_range[0]))
		except:
			notify_error.emit("There was an error in your Energy Levels window inputs")
			raise
		finally:
			self.plotting_finished.emit()
		

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
			self.fig.canvas.draw_idle()

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


peak_lock = threading.RLock()
def peak_locked(func):
	def _wrapper(*args, **kwargs):
		with peak_lock:
			return(func(*args, **kwargs))
	return _wrapper

class PeakfinderWindow(EQDockWidget):
	default_visible = False
	default_position = None
	available_in = ['LLWP',]

	peakfinding_started = pyqtSignal()
	peakfinding_finished = pyqtSignal()
	
	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		self.setWindowTitle("Peakfinder")

		layout = QVBoxLayout(margin=True)
		mainwidget = QWidget()
		mainwidget.setLayout(layout)
		self.setWidget(mainwidget)

		self.peaks = []

		class CustomPlotWidget(PlotWidget):
			def gui(self):
				super().gui()
				self.peaks_artist = self.ax.scatter([], [], color=config["peakfinder_peakcolor"], marker="*")

			def update_plot(self):
				xmin, xmax = self.xrange
				tuples = list(filter(lambda x: xmin < x[0] < xmax, self.parent.peaks))
				tuples = tuples if len(tuples)!=0 else [[None,None]]
				self.peaks_artist.set_offsets(tuples)
				self.peaks_artist.set_color(config["peakfinder_peakcolor"])
				super().update_plot()


		self.plot_widget = CustomPlotWidget(parent=self)
		layout.addWidget(self.plot_widget)

		tmp_layout = QHBoxLayout()
		layout.addLayout(tmp_layout)

		self.run_button = QQ(QPushButton, text="Run Peakfinder", change=lambda: self.find_peaks())
		tmp_layout.addWidget(self.run_button)
		tmp_layout.addWidget(QQ(QPushButton, text="Export Peaks", change=lambda: self.export_peaks()))

		uncert_input = QQ(QDoubleSpinBox, "peakfinder_width", range=(0, None), enabled=config["peakfinder_onlyunassigned"], minWidth=80)
		tmp_layout.addWidget(QQ(QCheckBox, "peakfinder_onlyunassigned", text="Only unassigned Lines", change=uncert_input.setEnabled))
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

			if key in config["peakfinder_kwargs"]:
				value = config["peakfinder_kwargs"][key]
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

		self.peakfinding_started.connect(lambda: self.run_button.setEnabled(False))
		self.peakfinding_started.connect(lambda: self.infolabel.setText("Finding peaks ..."))
		self.peakfinding_finished.connect(self.after_find_peaks)

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
		config["peakfinder_kwargs"] = kwargs
		return(kwargs.copy())

	@QThread.threaded_d
	@peak_locked
	def find_peaks(self, thread=None):
		self.peakfinding_started.emit()
		kwargs = self.get_kwargs()

		if "frequency" in kwargs:
			val = kwargs["frequency"]
			if type(val) in [tuple, list]:
				exp_df = ExpFile.get_data(xrange=val)
			else:
				exp_df = ExpFile.get_data().query(f"{val} < x")
			del kwargs["frequency"]
		else:
			exp_df = ExpFile.get_data()
		xs, ys = exp_df["x"].to_numpy(), exp_df["y"].to_numpy()

		peaks, props = signal.find_peaks(ys, **kwargs)
		self.peaks = np.array((xs[peaks], ys[peaks])).T
		self.plot_widget.update_plot()

		if config["peakfinder_onlyunassigned"]:
			assigned_xs = LinFile.get_data()["x"]
			uncertainty = config["peakfinder_width"]

			peaks_xs = self.peaks[:, 0]
			peaks_assigned = np.zeros(self.peaks.shape[0])

			for x in assigned_xs:
				peaks_assigned += (abs(peaks_xs - x) < uncertainty)
			self.peaks = self.peaks[peaks_assigned == 0]

		self.peaks = self.peaks[self.peaks[:, 1].argsort()[::-1]]
		self.peakfinding_finished.emit()

	@peak_locked
	def after_find_peaks(self):
		self.table.setRowCount(0)
		self.table.setColumnCount(2)
		self.table.setHorizontalHeaderLabels(["x", "y"])
		for x, y in self.peaks[:config["peakfinder_maxentries"]]:
			currRowCount = self.table.rowCount()
			self.table.insertRow(currRowCount)
			self.table.setItem(currRowCount, 0, QTableWidgetItem(f'{x:{config["flag_xformatfloat"]}}'))
			self.table.setItem(currRowCount, 1, QTableWidgetItem(f'{y:{config["flag_xformatfloat"]}}'))
		self.table.resizeColumnsToContents()
		self.table.setHidden(False)

		self.run_button.setEnabled(True)

		if len(self.peaks) > config["peakfinder_maxentries"]:
			self.infolabel.setText(f"Found {len(self.peaks)} peaks. Only the highest {config['peakfinder_maxentries']} entries are shown in the table for better performance. You can see all peaks by exporting or increasing the maximum number of displayed peaks in the configuration.")
		else:
			self.infolabel.setText(f"Found {len(self.peaks)} peaks.")

	@peak_locked
	def export_peaks(self):
		fname = QFileDialog.getSaveFileName(None, 'Choose file to save peaks to',"","CSV Files (*.csv);;All Files (*)")[0]
		if fname:
			np.savetxt(fname, self.peaks, delimiter="\t")

	@peak_locked
	def go_to(self, event):
		for idx in self.table.selectionModel().selectedIndexes():
			row_number = idx.row()

		if type(row_number) == int and len(self.peaks) > row_number:
			xmin, xmax = self.plot_widget.xrange
			width = (xmax - xmin)
			xcenter = self.peaks[row_number, 0]
			self.plot_widget.xrange = xcenter - width / 2 , xcenter + width / 2
			self.plot_widget.update_plot()

	def activateWindow(self):
		if hasattr(self, 'plot_widget'):
			self.plot_widget.from_current_plot()
		super().activateWindow()


cmd_lock = threading.RLock()
def cmd_locked(func):
	def _wrapper(*args, **kwargs):
		with cmd_lock:
			return(func(*args, **kwargs))
	return _wrapper

class CmdWindow(EQDockWidget):
	default_visible = False
	default_position = None
	available_in = ['LLWP', 'LASAP']

	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		self.setWindowTitle("Cmd Window")

		layout = QVBoxLayout(margin=True)
		mainwidget = QWidget()
		mainwidget.setLayout(layout)
		self.setWidget(mainwidget)

		self.tabs = QTabWidget()
		self.tabs.setTabsClosable(True)
		self.tabs.setMovable(True)
		self.tabs.setDocumentMode(True)

		if not config["cmd_commands"]:
			config["cmd_commands"] = [["Initial", "", True, True, True]]

		for values in config["cmd_commands"]:
			self.add_tab(*values)

		self.tabs.tabCloseRequested.connect(self.close_tab)
		self.tabs.tabBarDoubleClicked.connect(self.renameoradd_tab)
		self.tabs.setCurrentIndex(config["cmd_current"])
		self.tabs.currentChanged.connect(lambda x: self.update_cmd_command())

		layout.addWidget(QQ(QLabel, wordwrap=True, text="Here you can set up a command to run fitting and prediction programs to update your files."))
		layout.addWidget(self.tabs)
		layout.addWidget(QQ(QPushButton, text="Run", change=lambda x: self.run_pipe()))

	def add_tab(self, title, command=None, exp_checked=True, cat_checked=True, lin_checked=True):
		tmp = QWidget()
		layout = QVBoxLayout()

		if not command:
			command = ""

		layout.addWidget(QQ(QPlainTextEdit, value=command, change=lambda: self.update_cmd_command(), placeholder="Write your command line command here"))
		layout.addWidget(QQ(QCheckBox, text="Reread Exp Files after command finished", value=exp_checked, change=lambda x: self.update_cmd_command()))
		layout.addWidget(QQ(QCheckBox, text="Reread Cat Files after command finished", value=cat_checked, change=lambda x: self.update_cmd_command()))
		layout.addWidget(QQ(QCheckBox, text="Reread Lin Files after command finished", value=lin_checked, change=lambda x: self.update_cmd_command()))

		tmp.setLayout(layout)
		self.tabs.addTab(tmp, title)

	@cmd_locked
	def close_tab(self, index):
		tab = self.tabs.widget(index)
		tab.deleteLater()
		self.tabs.removeTab(index)
		if len(config["cmd_commands"]) > index:
			config["cmd_commands"].pop(index)
		if self.tabs.count() == 0:
			self.add_tab("New Tab")
			config["cmd_commands"].append(["New Tab", "", True, True, True])

	@cmd_locked
	def renameoradd_tab(self, index):
		if index == -1:
			self.add_tab("New Tab")
			config["cmd_commands"].append(["New Tab", "", True, True, True])
		elif self.tabs.widget(index) != 0:
			text, ok = QInputDialog().getText(self, "Tab Name","Enter the Tabs Name:")
			if ok and text:
				self.tabs.setTabText(index, text)
				config["cmd_commands"][index][0] = text

	@cmd_locked
	def update_cmd_command(self):
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

		config["cmd_commands"] = result
		config["cmd_current"] = self.tabs.currentIndex()

	@staticmethod
	@cmd_locked
	def run_pipe(index=None):
		if index == None:
			index = config["cmd_current"]

		if len(config["cmd_commands"]) == 0:
			notify_warning.emit("No Pipe command specified, therefore no Pipe process was started.")
			return

		title, command, exp_rr, cat_rr, lin_rr = config["cmd_commands"][index]

		if command == None:
			notify_warning.emit("No Pipe command specified, therefore no Pipe process was started.")
			return

		command = command.replace("\n", " && ")
		try:
			output = subprocess.check_output(command, shell=True)
			output = output.decode("utf-8")

			notify_info.emit(f"The subprocess was started and returned:\n{output}")

			for cls, reread in zip((ExpFile, CatFile, LinFile), (exp_rr, cat_rr, lin_rr)):
				if reread:
					cls.reread_all()

		except Exception as E:
			notify_error.emit(f"The command '{command}' failed with the Exception '{E}'.")
			raise

##
## Global Functions
##
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

	for _ in range(0, int(derivative)):
		ys = np.gradient(ys, edge_order=2)
	if derivative%2 == 0 and derivative != 0 and derivative%4 != 0:
		ys = -ys
	ymax = np.max(ys) if np.isfinite(ys).any() else 1
	if not np.isfinite(ymax) or ymax == 0:
		ymax = 1
	ys = amp*ys/ymax
	return(ys)

def fit_pgopher(xs, ys, peakdirection, fit_xs):
	ymin, ymax = np.min(ys), np.max(ys)

	if peakdirection < 0:
		cutoff = ymin + (ymax-ymin)/2
		mask = (ys <= cutoff)
	else:
		cutoff = ymax - (ymax-ymin)/2
		mask = (ys >= cutoff)

	fit_xs = xs[mask]
	fit_ys = ys[mask] - ymin

	xmiddle = np.sum(fit_xs*fit_ys)/np.sum(fit_ys)
	xuncert = 0

	return(xmiddle, xuncert, fit_xs, fit_ys + ymin)

def fit_polynom(xs, ys, peakdirection, fit_xs, rank):
	try:
		popt = np.polyfit(xs, ys, rank)
	except Exception as E:
		popt = np.polyfit(xs, ys, rank)
	polynom = np.poly1d(popt)
	fit_ys = polynom(fit_xs)

	if peakdirection < 0:
		xmiddle = fit_xs[np.argmin(fit_ys)]
	else:
		xmiddle = fit_xs[np.argmax(fit_ys)]

	xuncert = 0
	return(xmiddle, xuncert, fit_xs, fit_ys)

def fit_polynom_multirank(xs, ys, peakdirection, fit_xs, maxrank):
	best_rms = np.inf
	best_rank = 0

	maxrank = min(len(xs), maxrank)
	for rank in range(maxrank):
		try:
			try:
				popt = np.polyfit(xs, ys, rank)
			except Exception as E:
				popt = np.polyfit(xs, ys, rank)
			polynom = np.poly1d(popt)
			fit_ys = polynom(exp_xs)

			rms = np.mean((fit_ys - ys)**2)
			if rms < best_rms:
				best_rms = rms
				best_rank = rank
		except Exception as E:
			continue

	popt = np.polyfit(xs, ys, best_rank)
	polynom = np.poly1d(popt)
	fit_ys = polynom(fit_xs)

	if peakdirection < 0:
		xmiddle = fit_xs[np.argmin(fit_ys)]
	else:
		xmiddle = fit_xs[np.argmax(fit_ys)]

	xuncert = 0
	return(xmiddle, xuncert, fit_xs, fit_ys)

def fit_lineshape(xs, ys, peakdirection, fit_xs, profilname, derivative, offset, **kwargs):
	xmin, xmax = xs.min(), xs.max()
	x0 = (xmin + xmax) / 2
	
	xs_weight_factor = kwargs.get('xs_weight_factor', 4)
	if xs_weight_factor:
		ys_weighted = ys * np.exp(- np.abs(np.abs(xs - x0) / (xmax - xmin) ) * xs_weight_factor)
		x0 = xs[np.argmax(ys_weighted)] if peakdirection >= 0 else xs[np.argmin(ys_weighted)]
	
	ymin, ymax, ymean, yptp = ys.min(), ys.max(), ys.mean(), np.ptp(ys)
	y0 = 0
	
	w0 = kwargs.get('w0', (xmax - xmin) / 10 )
	wmin = kwargs.get('wmin', 0)
	wmax = kwargs.get('wmax', (xmax - xmin) )
	
	amp_min, amp_max = -3*yptp, 3*yptp
	if peakdirection < 0:
		amp_max = 0
		y0 = -yptp
	if peakdirection > 0:
		amp_min = 0
		y0 = yptp

	p0 = [x0, y0, w0] if profilname != 'Voigt' else [x0, y0, w0, w0]
	bounds = [[xmin, amp_min, wmin], [xmax, amp_max, wmax]] if profilname != 'Voigt' else [[xmin, amp_min, wmin, wmin], [xmax, amp_max, wmax, wmax]]
	function = lambda *x: lineshape(profilname, derivative, *x)

	if offset:
		function = lambda *x: lineshape(profilname, derivative, *x[:-1]) + x[-1]
		p0.append(ymean)
		bounds[0].append(ymin)
		bounds[1].append(ymax)

	try:
		popt, pcov = optimize.curve_fit(function, xs, ys, p0=p0, bounds=bounds)
	except Exception as E:
		popt, pcov = optimize.curve_fit(function, xs, ys, p0=p0, bounds=bounds)
	perr = np.sqrt(np.diag(pcov))
	fit_ys = function(fit_xs, *popt)

	xmiddle = popt[0]
	xuncert = perr[0]

	return(xmiddle, xuncert, fit_xs, fit_ys)

def get_fitfunction(fitmethod, offset=False, **kwargs):
	fit_function = {
		'Pgopher': fit_pgopher,
		'Polynom': lambda *args: fit_polynom(*args, config['fit_polynomrank']),
		'MultiPolynom': lambda *args: fit_polynom_multirank(*args, config['fit_polynommaxrank']),
	}.get(fitmethod)
	
	if not fit_function:
		profilname, derivative = {
			'Gauss':					('Gauss', 0),
			'Lorentz':					('Lorentz', 0),
			'Voigt':					('Voigt', 0),
			'Gauss 1st Derivative':		('Gauss', 1),
			'Lorentz 1st Derivative':	('Lorentz', 1),
			'Voigt 1st Derivative':		('Voigt', 1),
			'Gauss 2nd Derivative':		('Gauss', 2),
			'Lorentz 2nd Derivative':	('Lorentz', 2),
			'Voigt 2nd Derivative':		('Voigt', 2),
		}[fitmethod]
		fit_function = lambda *args, kwargs=kwargs: fit_lineshape(*args, profilname, derivative, offset, **kwargs)
	return(fit_function)

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

def symmetric_ticklabels(ticks):
	float_format = config['flag_xformatfloat']
	tick_labels = []
	for a, o in zip(ticks, ticks[::-1]):
		if not (np.isfinite(a) and np.isfinite(o)):
			continue
		dec_a = len(f"{a:{float_format}}".rstrip("0").split(".")[1])
		dec_o = len(f"{o:{float_format}}".rstrip("0").split(".")[1])
		if dec_a == dec_o:
			tick_labels.append(f"{a:{float_format}}".rstrip("0").rstrip("."))
		else:
			trailing_zeros = 4 - max(dec_a, dec_o)
			tick = f"{a:{float_format}}"[:-trailing_zeros] if trailing_zeros else f"{a:{float_format}}"
			tick_labels.append(tick)
	return(tick_labels)

def is_dark_theme():
	return(QApplication.styleHints().colorScheme() == Qt.ColorScheme.Dark)

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
	if kwargs.get("buttons", True) is False:
		widget.setButtonSymbols(QAbstractSpinBox.ButtonSymbols.NoButtons)

	if widgetclass in [QSpinBox, QDoubleSpinBox, QDoubleSpinBoxFullPrec]:
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
		setter(config[config_key])
		changer(lambda x=None, key=config_key: config.__setitem__(key, getter(), widget))
		config.register_widget(config_key, widget, lambda: setter(config[config_key]))
	if "change" in kwargs:
		changer(kwargs["change"])
	if "changes" in kwargs:
		for change in kwargs["changes"]:
			changer(change)

	return widget

def bin_data(dataframe, binwidth, range):
	# It is not verified that the data to be binned is
	#  - from the same file (-> group by files)
	#  - is equidistant due to filters or different files being merged
	#    (-> bins can have different number of points)
	# Otherwise a much faster algorithm could be implemented

	## The old version (commented out) is significantly slower than the new version.
	## Tested with cyclopentadiene spectrum around 200 GHz

	# Width [MHz]| Old [ms] | New [ms] 
	# -----------|----------|----------
	#          1 |     1.73 |     1.42 
	# 		  10 |     1.07 |     1.08 
	# 	     100 |     1.14 |     1.01 
	# 	    1000 |     2.64 |     1.67 
	#      10000 |    18.35 |     7.38 
	#     100000 |   160.05 |    53.30 
	#    1000000 |   569.59 |   163.29 
	#   10000000 |   618.89 |   171.58 
	#  100000000 |   668.12 |   191.47 

	# length = len(dataframe)
	# dataframe.loc[:,"bin"] = (dataframe.loc[:,"x"]-range[0]) // binwidth
	# # For assignments (lin_df) as they do not have an intensity
	# if "y" not in dataframe:
	# 	dataframe = dataframe.loc[dataframe.drop_duplicates(("bin", "filename"), keep="last").sort_values(["x"]).index]
	# else:
	# 	dataframe = dataframe.loc[dataframe.sort_values("y").drop_duplicates(("bin", "filename"), keep="last").sort_values(["x"]).index]
	# return(dataframe)

	length = len(dataframe)
	dataframe.loc[:,"bin"] = (dataframe.loc[:,"x"]-range[0]) // binwidth
	
	if "y" not in dataframe:
		index_ = dataframe.groupby(['bin', 'filename'], observed=True)['x'].idxmax()
	else:
		index_ = dataframe.groupby(['bin', 'filename'], observed=True)['y'].idxmax()
	
	dataframe = dataframe.loc[index_]
	return(dataframe)


def llwpfile(extension=""):
	home = os.path.expanduser("~")
	llwpfolder = os.path.join(home, f".{APP_TAG.lower()}")
	
	if not os.path.isdir(llwpfolder):
		os.mkdir(llwpfolder)

	return(os.path.join(llwpfolder, extension))

def exp_to_df(fname, n_bytes=4096, **kwargs):
	kwargs = {
		'dtype': np.float64,
		'names': ['x', 'y'],
		'usecols': [0, 1],
		'comment': '#',
		'header': None,
		'engine': 'c',
		** kwargs,
	}

	# Find the delimiter automatically if not specified
	if 'delimiter' not in kwargs and 'sep' not in kwargs:
		sniffer = csv.Sniffer()
		with open(fname, 'r') as file:
			data = file.read(n_bytes)
		delimiter = sniffer.sniff(data).delimiter
		kwargs['delimiter'] = delimiter

	data = pd.read_csv(fname, **kwargs)
	data = data.dropna()
	data['filename'] = fname

	return(data)

def restart():
	project_filename = llwpfile('.files')
	File.save_files(project_filename)
	mainwindow.closeEvent()
	NewAssignments.get_instance().save_backup()
	os.execl(sys.executable, sys.executable, __file__, project_filename)

def except_hook(cls, exception, traceback):
	if isinstance(exception, GUIAbortedError):
		return
	
	if issubclass(cls, KeyboardInterrupt):
		sys.exit(0)

	sys.__excepthook__(cls, exception, traceback)
	with open(llwpfile(".err"), "a+", encoding="utf-8") as file:
		time_str = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
		file.write(f"{time_str}: \n{exception}\n{''.join(tb.format_tb(traceback))}\n\n")
	try:
		NewAssignments.get_instance().save_backup()
		notify_error.emit(f"{exception}\n{''.join(tb.format_tb(traceback))}")
	except Exception as E:
		pass

def get_all_subclasses(cls):
	all_subclasses = []
	for subclass in cls.__subclasses__():
		all_subclasses.append(subclass)
		all_subclasses.extend(get_all_subclasses(subclass))
	return all_subclasses




# ASAP
class SpecialFilesHandlerASAP():
	@classmethod
	def save_files(cls, dict_):
		if ASAPAx.egy_filename is not None:
			dict_['.egy'] = ASAPAx.egy_filename
		return(dict_)
	
	@classmethod
	def load_files(cls, dict_):
		if '.egy' in dict_:
			egy_filename = dict_.pop('.egy')
			ASAPAx.load_egy_file(egy_filename)
		return(dict_)
	
	@classmethod
	def sort_files_by_type(cls, files):
		egy_files = []
		other_files = []

		for file in files:
			is_egy_file = (os.path.splitext(file)[1] == '.egy')

			if is_egy_file:
				egy_files.append(file)
			else:
				other_files.append(file)
		
		if len(egy_files) > 1:
			notify_warning.emit('Only a single *.egy file can be loaded at the same time.')
		elif len(egy_files) == 1:
			ASAPAx.load_egy_file(egy_files[0])
		return(other_files)

class LevelSelector(SeriesSelector):
	values_changed = pyqtSignal()

	def __init__(self, parent, initial_values={}):
		super(SeriesSelector, self).__init__(parent)
		self.n_qns = N_QNS
		self.parent = parent
		self.updating = False
		self.state = {
			'is_upper_state': True,
			'qns': (1, 0, 1, 0, 0, 0),
			'incr': (True, False, True, False, False, False),
			'diff': (1, 0, 1, 0, 0, 0),
			'use_diff': False,
		}
		self.state.update(initial_values)

		layout = QGridLayout()

		create_qn = lambda: QQ(QSpinBox, minWidth=60, maxWidth=60, range=(None, None), visible=False,
							singlestep=1, change=lambda x: self.changed())

		self.qns = [create_qn() for _ in range(self.n_qns)]

		create_widget = lambda widget, kwargs: QQ(widget, minWidth=40, maxWidth=40, visible=False,
							change=lambda x: self.changed(), **kwargs)

		self.incr = [create_widget(QCheckBox, {'text':'Inc'}) for _ in range(self.n_qns)]
		self.diff = [create_widget(QSpinBox, {'range': (None, None), 'singlestep': 1, }) for _ in range(self.n_qns)]

		self.incqns = QQ(QPushButton, text="Inc", change=lambda x: self.incdecqns(+1), width=40)
		self.decqns = QQ(QPushButton, text="Dec", change=lambda x: self.incdecqns(-1), width=40)

		self.togglediff = QQ(QToolButton, text="⇆", change=lambda x: self.change_incr_mode(), width=40)
		
		for i, widget in enumerate(self.qns):
			layout.addWidget(widget, 0, i)

		for i, incr, diff in zip(range(N_QNS), self.incr, self.diff):
			tmp = QHBoxLayout()
			tmp.addWidget(incr)
			tmp.addWidget(diff)

			layout.addLayout(tmp, 4, i)
			# layout.setColumnStretch(i, 100)

		# for i in range(self.n_qns):
		# 	# layout.setColumnStretch(i, 100)
		
		layout.addWidget(self.togglediff, 5, self.n_qns, 1, 1)

		layout.addWidget(self.incqns, 0, self.n_qns, 1, 2)
		layout.addWidget(self.decqns, 4, self.n_qns, 1, 2)
		
		self.is_upper_state_checkbox = QQ(QCheckBox, change=lambda x: self.changed(), text='Upper State Level')
		layout.addWidget(self.is_upper_state_checkbox, 5, 0, 1, self.n_qns)

		layout.setRowStretch(6, 10)
		layout.setColumnStretch(self.n_qns+2, 1)

		self.layout = layout
		self.setLayout(layout)

		self.set_state()

		config.register_widget("series_qns", self.togglediff, self.set_state)
		config.register_widget("flag_showseriesarrows", self.togglediff, self.change_arrows)
		self.change_arrows()

	def set_state(self):
		self.updating = True
		state = self.state

		n_qns = config["series_qns"]
		are_visible = [True] * n_qns +[False] * (self.n_qns - n_qns)

		for qn, value, is_visible in zip(self.qns, state["qns"], are_visible):
			qn.setValue(value)
			qn.setVisible(is_visible)
		for widget, value, is_visible in zip(self.incr, state["incr"], are_visible):
			widget.setChecked(value)
			widget.setVisible(is_visible and not state['use_diff'])
		for widget, value, is_visible in zip(self.diff, state["diff"], are_visible):
			widget.setValue(value)
			widget.setVisible(is_visible and state['use_diff'])
		
		self.is_upper_state_checkbox.setChecked(state['is_upper_state'])
		self.updating = False
		self.changed()

	def incdecqns(self, dir):
		self.updating = True
		incr_values = (x.value() for x in self.diff) if self.state['use_diff'] else (x.isChecked() for x in self.incr)

		for qn, incr in zip(self.qns, incr_values):
			qn.setValue(qn.value()+dir*incr)
		self.updating = False
		self.changed()

	def changed(self):
		if self.updating:
			return
	
		self.state["qns"] = [x.value() for x in self.qns]
		self.state["incr"] = [x.isChecked() for x in self.incr]
		self.state["diff"] = [x.value() for x in self.diff]
		self.state["is_upper_state"] = self.is_upper_state_checkbox.isChecked()
		self.values_changed.emit()
	
	def change_arrows(self):
		show_arrows = config['flag_showseriesarrows']
		width = 60 if show_arrows else 40
		button_symbols = QAbstractSpinBox.ButtonSymbols.UpDownArrows if show_arrows else QAbstractSpinBox.ButtonSymbols.NoButtons

		for widget in self.qns + self.diff:
			widget.setMinimumWidth(width)
			widget.setMaximumWidth(width)
			widget.setButtonSymbols(button_symbols)

class ASAPAx(LWPAx):
	fit_vline = None
	fit_curve = None
	fit_methods = ('Pgopher', 'Polynom', 'MultiPolynom', 'Gauss', 'Lorentz')

	egy_df = None
	egy_filename = None

	def __init__(self, ax, row_i, col_i):
		self.ax = ax
		self.row_i = row_i
		self.col_i = col_i
		
		self.xrange = (-1, 1)
		self.annotation = None
		self.qns = None
		self.entries = None
		self.is_upper_state = True

		self.corr_xs = None
		self.corr_ys = None
		self.curr_state = {}
		
		with matplotlib_lock:
			self.span =  matplotlib.widgets.SpanSelector(ax, lambda xmin, xmax: self.on_range(xmin, xmax), 'horizontal', useblit=True, button=1)
		
		self.exp_coll = matplotlib.collections.LineCollection(np.zeros(shape=(0,2,2)), colors=config["color_exp"], capstyle='round')

		with matplotlib_lock:
			ax.add_collection(self.exp_coll)
			
			ax.yaxis.set_visible(False)
			if row_i:
				ax.xaxis.set_visible(False)
			else:
				ax.set_xticks([])

			if row_i != 0:
				ax.spines['bottom'].set_visible(False)
			if row_i != config['plot_rows'] - 1:
				ax.spines['top'].set_visible(False)
			if col_i != 0:
				ax.spines['left'].set_visible(False)
			if col_i != config['plot_cols'] - 1:
				ax.spines['right'].set_visible(False)

	def vals_to_coll(self, xs, ys):
		if xs is None or ys is None:
			segs = np.array([])
		else:
			xrange = self.xrange
			bins = config['plot_bins']
			nobinning = config['plot_skipbinning']
			binwidth = (xrange[1]-xrange[0]) / bins

			if len(xs) > max(bins, nobinning) and binwidth != 0:
				df = pd.DataFrame({'x': xs, 'y': ys, 'filename': 0})
				df = bin_data(df, binwidth, xrange)
				xs, ys = df['x'], df['y']

			segs = np.array(((xs[:-1], xs[1:]), (ys[:-1], ys[1:]))).T

		self.exp_coll.set(segments=segs, color=config['color_exp'])

	# @status_d
	@QThread.threaded_d
	@drawplot_decorator.d
	def update(self, thread=None):
		ax = self.ax

		offset, width = config['plot_offset'], config['plot_width']
		self.xrange = (offset - width/2, offset + width/2)
		
		ax.set_xlim(self.xrange)
		self.set_xticklabels()
		yrange = [-1, 1]

		tot_xs = self.corr_xs
		tot_ys = self.corr_ys


		if tot_xs is not None and len(tot_xs):
			min_index = tot_xs.searchsorted(self.xrange[0], side='right')
			max_index = tot_xs.searchsorted(self.xrange[1], side='left')

			tot_xs = tot_xs[min_index:max_index]
			tot_ys = tot_ys[min_index:max_index]
		
		
		self.vals_to_coll(tot_xs, tot_ys)

		if tot_ys is not None and len(tot_ys):
			yrange = (np.min(tot_ys), np.max(tot_ys))
		margin = config['plot_ymargin']

		

		yrange = (yrange[0]-margin*(yrange[1]-yrange[0]), yrange[1]+margin*(yrange[1]-yrange[0]))
		if np.isnan(yrange[0]) or np.isnan(yrange[1]) or yrange[0] == yrange[1]:
			yrange = (-2,+2)

		ax.set_ylim(yrange)
		
		self.update_annotation()
		

	def update_annotation(self):
		fstring = config['plot_annotationfstring']

		if not fstring:
			if self.annotation:
				self.annotation.remove()
				self.annotation.set_visible(False)
				self.annotation = None
			return

		color = matplotlib.rcParams['text.color']

		if self.qns is not None:
			qnstring = ','.join(map(str, self.qns))

			query = ' and '.join([f'(qnu{i+1} == {qn} and qnl{i+1} == 0)' for i, qn in enumerate(self.qns)])
			if query and LinFile.has_results(query):
				color = config["color_lin"]
		else:
			qnstring = ''


		vars = {
			'x': config['plot_offset'],
			'qns': qnstring,
			'width': config['plot_width'],
		}
		text = fstring.format(**vars)

		if self.annotation is None:
			kwargs = config['plot_annotationkwargs']
			ax = self.ax
			self.annotation = ax.text(**kwargs, s=text, color=color, transform=ax.transAxes)
		else:
			if self.annotation.get_text() != text:
				self.annotation.set_text(text)
			if matplotlib.colors.to_rgba(self.annotation.get_color()) != matplotlib.colors.to_rgba(color):
				self.annotation.set_color(color)

	def set_xticklabels(self):
		if self.row_i:
			return
		
		ax = self.ax
		ticks = np.linspace(*self.xrange, config['plot_xticks'])
		tickformat = config['plot_xtickformat']

		if tickformat == 'scientific':
			ticklabels = [f"{x:.2e}".replace("e+00", "").rstrip("0").rstrip(".") for x in ticks]
		else:
			ticklabels = symmetric_ticklabels(ticks)		

		if self.col_i and len(ticks) > 1:
			ticks = ticks[1:]
			ticklabels = ticklabels[1:]
		ax.set_xticks(ticks)
		ax.set_xticklabels(ticklabels)
	
	def create_qns_dict(self, complete=False):
		qns_dict = {} if not complete else {f'qn{ul}{i+1}': pyckett.SENTINEL for ul in 'ul' for i in range(N_QNS)}
		if self.qns is None:
			return(qns_dict)
		
		for i, qn in enumerate(self.qns):
			qns_dict[f'qnu{i+1}'] = qn
			qns_dict[f'qnl{i+1}'] = 0
		return(qns_dict)

	def on_range(self, xmin, xmax):
		xmin_ax, xmax_ax = self.xrange
		if xmax == xmin or xmax > xmax_ax or xmin < xmin_ax:
			return
		self.fit_data(xmin, xmax)

	def fit_data(self, xmin, xmax, onclick=False):
		if self.qns is None:
			return

		# Delete artists highlighting previous fit
		if self.__class__.fit_vline is not None:
			self.__class__.fit_vline.remove()
			self.__class__.fit_vline = None
		if self.__class__.fit_curve is not None:
			self.__class__.fit_curve.remove()
			self.__class__.fit_curve = None
		
		if onclick:
			xmiddle = xmin
			xuncert = 0
		
		else:
			# Fit the data
			xmiddle, xuncert, fit_xs, fit_ys = self.fit_peak(xmin, xmax)
			self.__class__.fit_curve = self.ax.plot(fit_xs, fit_ys, color=config["color_fit"], alpha=0.7, linewidth=1)[0]
		
		# Highlight fit in plot
		self.__class__.fit_vline = self.ax.axvline(x=xmiddle, color=config["color_fit"], ls="--", alpha=1, linewidth=1)

		# Emit signal for detail viewer
		mainwindow.lwpwidget.peak_fitted.emit(self, xmiddle)

		# Add predicted energy to offset
		qns_dict = self.create_qns_dict(complete=True)
		egy_df = self.__class__.egy_df
		egy_val = 0
		if egy_df is not None:
			query = ' and '.join([f'qn{i+1} == {qn}' for i, qn in enumerate(self.qns)])
			vals = egy_df.query(query)["egy"].to_numpy()
			egy_val = vals[0] if len(vals) else 0
		
		if egy_val == 0:
			notify_warning.emit('No corresponding energy level found! Please check if an energy file is loaded.')

		error = self.fit_determine_uncert(0, xmiddle, xuncert)
		
		if self.is_upper_state:
			xenergy = xmiddle + egy_val
		else:
			xenergy = xmiddle - egy_val

		# Create assignment object
		new_assignment = {'x': xenergy, 'error': error, 'xpre': 0}
		new_assignment.update(qns_dict)
		new_assignment.update({'weight': 1, 'comment': config['fit_comment'], 'filename': '__newassignments__'})
		
		if not config['asap_assigntransitions']:
			NewAssignments.get_instance().add_row(new_assignment)
		
		else:
			# Entries are already in energy units -> same error as energy level is appropriate
			entries = self.entries.query('(use_for_cross_correlation)').copy()
			entries['x'] += xmiddle
			noq = config['series_qns']
			energy_assignment_df = pd.DataFrame(new_assignment, index=[0])
			entries = pd.concat( (energy_assignment_df, entries) )
			
			new_assignments = {
				'xpre': 0,
				'x': entries['x'],
				'weight': 1,
				'error': error,
				'comment': config['fit_comment'],
				'filename': '__newassignments__',
			}

			for i in range(noq):
				new_assignments[f'qnu{i+1}'] = entries[f'qnu{i+1}']
				new_assignments[f'qnl{i+1}'] = entries[f'qnl{i+1}']
			
			for i in range(noq, N_QNS):
				new_assignments[f'qnu{i+1}'] = pyckett.SENTINEL
				new_assignments[f'qnl{i+1}'] = pyckett.SENTINEL

			NewAssignments.get_instance().add_rows(new_assignments)
		

	def fit_peak(self, xmin, xmax):
		exp_xs, exp_ys = self.corr_xs, self.corr_ys
		min_index = exp_xs.searchsorted(xmin, side='right')
		max_index = exp_xs.searchsorted(xmax, side='left')

		exp_xs, exp_ys = exp_xs[min_index:max_index], exp_ys[min_index:max_index]
		
		fit_xs = np.linspace(xmin, xmax, config['fit_xpoints'])

		peakdirection = config['fit_peakdirection']
		fitmethod = config['fit_fitmethod']

		if (len(exp_xs) == 0) or ((len(exp_xs) < 2) and fitmethod != 'Pgopher'):
			notify_error.emit('The data could not be fit as there were too few points selected.')
			raise GUIAbortedError('The data could not be fit as there were too few points selected.')

		try:
			fit_function = get_fitfunction(fitmethod, config['fit_offset'])
			xmiddle, xuncert, fit_xs, fit_ys = fit_function(exp_xs, exp_ys, peakdirection, fit_xs)

		except Exception as E:
			self.fitcurve = None
			self.fitline = None
			notify_error.emit(f"The fitting failed with the following error message : {str(E)}")
			raise

		return(xmiddle, xuncert, fit_xs, fit_ys)


	@classmethod
	def load_egy_file(cls, filename=None):
		if not filename:
			filename, filter = QFileDialog.getOpenFileName(None, 'Choose *.egy file')
			if not filename:
				return
		
		egy_df = pyckett.egy_to_df(filename, sort=False)
		cls.egy_df = egy_df
		cls.egy_filename = filename
		
		basename = os.path.basename(filename)

		if hasattr(ASAPSettingsWindow, 'instance'):
			ASAPSettingsWindow.instance.egy_file_button.setText(basename)
		
		notify_info.emit(f'Successfully loaded the energy file \'{basename}\'.')

	@classmethod
	def fit_determine_uncert(cls, *args, **kwargs):
		error = super().fit_determine_uncert(*args, **kwargs)
		return(-abs(error))

		
class ASAPMenu(Menu):
	def __init__(self, parent, *args, **kwargs):
		mb = parent.menuBar()
		
		# Create top level menus
		top_menu_labels = ("Files", "View", "Fit", "Info")
		self.top_menus = {}
		for label in top_menu_labels:
			menu = mb.addMenu(f'{label}')
			self.top_menus[label] = menu
		

		toggleaction_files = FileWindow.instance.toggleViewAction()
		toggleaction_files.setText('Edit Files')
		toggleaction_files.setShortcut('Shift+1')

		toggleaction_config = ConfigWindow.instance.toggleViewAction()
		toggleaction_config.setShortcut('Shift+0')

		toggleaction_credits = CreditsWindow.instance.toggleViewAction()
		toggleaction_credits.setText("Credits and License")
		toggleaction_credits.setToolTip("See the Credits and License")


		view_actions = [
			ASAPSettingsWindow.instance.toggleViewAction(),
			NewAssignmentsWindow.instance.toggleViewAction(),
			ASAPDetailViewer.instance.toggleViewAction(),
			LogWindow.instance.toggleViewAction(),
			CmdWindow.instance.toggleViewAction(),
		]
		
		for i, view_action in enumerate(view_actions):
			view_action.setShortcut(f'Shift+{i+2}')


		fitfunction_menu = QMenu("Choose Fit Function", parent=parent)
		self.fitfunction_actions = {}

		current_method = config['fit_fitmethod']
		for method in ASAPAx.fit_methods:
			is_checked = (method == current_method)
			callback = lambda _, method=method: self.set_fitmethod(method)
			self.fitfunction_actions[method] = QQ(QAction, parent=parent, text=f"{method}", change=callback, checkable=True, value=is_checked)
			fitfunction_menu.addAction(self.fitfunction_actions[method])
		config.register('fit_fitmethod', self.on_fitfunction_changed)

		actions = {
			'Files': (
				QQ(QAction, parent=parent, text="Add Files", change=File.add_files_dialog, shortcut='Ctrl+O', tooltip="Add any kind of Files"),
				QQ(QAction, parent=parent, text='Reread All Files', change=lambda _: File.reread_all(), shortcut='Ctrl+R', tooltip="Reread all Exp, Cat and Lin files"),
				None,
				toggleaction_files,
				None,
				QQ(QAction, parent=parent, text="Save current values as default", tooltip="Save current configuration as default", change=lambda _: config.save()),
				None,
				QQ(QAction, parent=parent, text="Save Files as Project", change=lambda _: File.save_files_gui(), tooltip="Save all loaded files and their parameters as a project."),
				QQ(QAction, parent=parent, text="Load Project", change=lambda _: File.load_files_gui(), tooltip="Load a project."),
				None,
				QQ(QAction, parent=parent, text="Quit", change=mainwindow.close),
			),
			'Fit': (
				fitfunction_menu,
				QQ(QAction, parent=parent, text="Change Function", shortcut="Ctrl+F", tooltip="Cycle through the available fit-functions", change=lambda _: self.next_fitmethod()),
				None,
				QQ(QAction, parent=parent, text="Change Fit Color", tooltip="Change the color of the fitfunction", change=lambda _: self.change_fitcolor()),
			),
			'View': (
				toggleaction_config,
				None,
				*view_actions,
				None,
				QQ(QAction, parent=parent, text='Open Console', shortcut='CTRL+K', change= lambda _: ConsoleDialog.show_dialog()),
				None,
			),
			'Info': (
				QQ(QAction, parent=parent, text="Open LASAP folder", change=lambda x: webbrowser.open(f'file:///{llwpfile()}'), tooltip="Open the folder containing the config, ...", ),
				QQ(QAction, parent=parent, text="Send Mail to Author", tooltip="Send a mail to the developer", change=lambda x: self.send_mail_to_author()),
				toggleaction_credits,
			)
			
		}

		for label, menu in self.top_menus.items():
			for widget in actions.get(label, []):
				if widget is None:
					menu.addSeparator()
				elif isinstance(widget, QAction):
					menu.addAction(widget)
				else:
					menu.addMenu(widget)
		
	def next_fitmethod(self):
		fitmethods = ASAPAx.fit_methods
		newindex = fitmethods.index( config['fit_fitmethod'] ) + 1
		newindex = newindex % len(fitmethods)
		newvalue = fitmethods[newindex]
		config['fit_fitmethod'] = newvalue

class ASAPWidget(LWPWidget):
	_ax_class = ASAPAx
	peak_fitted = pyqtSignal(ASAPAx, float)


	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		
		elem = QQ(QLabel, text='    Offset: ', visible=config['isvisible_controlswidth'])
		config.register("isvisible_controlswidth", lambda elem=elem: elem.setVisible(config["isvisible_controlswidth"]))
		self.toplayout.addWidget(elem)
		
		elem = QQ(QDoubleSpinBox, 'plot_offset', range=(None, None), minWidth=85, visible=config['isvisible_controlswidth'])
		config.register("isvisible_controlswidth", lambda elem=elem: elem.setVisible(config["isvisible_controlswidth"]))
		self.toplayout.addWidget(elem)


	def prefilter_correlation_plot_entries(self, state, n_rows, n_qns, offset, width):
		qns = state['qns'][:n_qns]
		diffs = (state['diff'] if state['use_diff'] else state['incr'])[:n_qns]

		qns = np.array(qns)
		diffs = np.array(diffs)

		# Prefilter df to all transitins belonging to series
		conditions, conditions_incr = ['(visible)'], []
		normalizing_value = None
		ul = 'u' if state['is_upper_state'] else 'l'
		
		for i, qn, diff in zip(range(n_qns), qns, diffs):
			diff = int(diff)
			if diff:
				if normalizing_value is None:
					normalizing_value = qn // diff
				conditions_incr.append(f"((qn{ul}{i+1} - {qn-normalizing_value*diff})/{diff})")
			else:
				conditions.append(f"(qn{ul}{i+1} == {qn})")

		if len(conditions_incr):
			conditions.append(" == ".join(conditions_incr))
		
		conditions = " and ".join(conditions)
		entries = CatFile.query_c(conditions).copy()
		
		if config['asap_catunitconversionfactor']:
			entries['x'] *= config['asap_catunitconversionfactor']

		if config['asap_query']:
			entries['use_for_cross_correlation'] = entries.eval(config['asap_query'])
		else:
			entries['use_for_cross_correlation'] = True

		entries['xmin'] = entries['x'] + offset - width/2
		entries['xmax'] = entries['xmin'] + width
		entries['min_index'], entries['max_index'] = ExpFile.xs_to_indices(entries['xmin'], entries['xmax'])

		return(entries, qns, diffs, ul)

	def calc_correlation_plot(self, row_qns, row_entries, offset, width, resolution):
		transitions = row_entries[row_entries['use_for_cross_correlation']]

		ref_xs = transitions['x']
		ref_ys = transitions['y']

		min_indices, max_indices = transitions['min_index'], transitions['max_index']

		xmin, xmax = offset - width/2, offset + width/2
		tot_xs = np.arange(xmin, xmax + resolution, resolution)
		tot_ys = np.ones_like(tot_xs)

		exp_len = len(ExpFile.df)
		n_correlated_transitions = 0
		use_weights = config['asap_weighted']

		minimum_intensity = np.log10(ref_ys.min())

		for min_index, max_index, ref_pos, ref_int in zip(min_indices, max_indices, ref_xs, ref_ys):
			min_index = max(0, min_index-1)
			max_index = min(exp_len, max_index+1)

			dataframe = ExpFile.df.iloc[min_index:max_index].copy()
			dataframe = dataframe[dataframe['visible']]

			if not len(dataframe):
				continue
			
			xs = dataframe['x']
			ys = dataframe['y']

			interp_ys = np.interp(tot_xs, xs-ref_pos, ys, left=1, right=1)

			# @Luis: Check this
			# Christian uses the minimal intensity from query as the minimum_intensity
			if use_weights:
				power = np.log10(ref_int) - minimum_intensity + 1
				interp_ys = np.power(np.abs(interp_ys), power)
				
			tot_ys *= interp_ys
			n_correlated_transitions += 1

		if n_correlated_transitions < 2:
			tot_xs = tot_ys = np.array([])
			
		# @Luis: Think about offering options to normalize spectrum here

		return(tot_xs, tot_ys)


	@QThread.threaded_d
	@status_d
	@drawplot_decorator.d
	def calc_correlation_plots(self, thread=None):
		if not hasattr(ASAPSettingsWindow, 'instance'):
			return
		if not hasattr(ASAPSettingsWindow.instance, 'tab'):
			return

		n_rows = config['plot_rows']
		n_cols = config['plot_cols']
		n_qns = config['series_qns']

		thread.earlyreturn()

		# Calculate positions and qns for each column
		tab_widget = ASAPSettingsWindow.instance.tab
		n_widgets = tab_widget.count()

		offset = config['plot_offset']
		width = config['plot_width']
		resolution = config['asap_resolution']

		threads = []
		with matplotlib_lock:
			if self.lwpaxes.shape != (n_rows, n_cols):
				notify_error.emit('Shape of LWPAxes is out of sync with requested values.')
				return

			for i_col in range(n_widgets):
				thread.earlyreturn()
				if i_col > n_cols:
					continue

				refwidget = tab_widget.widget(i_col)
				state = refwidget.state
				
				entries, qns, diffs, ul = self.prefilter_correlation_plot_entries(state, n_rows, n_qns, offset, width)
				
				thread.earlyreturn()

				# Get the specific entries for each ax 
				for i_row in range(n_rows):
					row_qns = qns + i_row * diffs
					cond = [f"(qn{ul}{i+1} == {qn})" for i, qn in enumerate(row_qns)]
					condition  = " & ".join(cond)
					row_entries_all = entries.query(condition)
					row_entries = row_entries_all[row_entries_all['use_for_cross_correlation']]

					corr_xs, corr_ys = self.calc_correlation_plot(row_qns, row_entries, offset, width, resolution)
					
					ax = self.lwpaxes[i_row, i_col]
					ax.entries = row_entries_all
					ax.corr_xs = corr_xs
					ax.corr_ys = corr_ys
					ax.qns = row_qns
					ax.is_upper_state = state['is_upper_state']
					threads.append(ax.update())
					
					thread.earlyreturn()

			# Edge cases of no reference tabs and too few tabs
			for i_col in range(n_widgets, n_cols):
				for i_row in range(n_rows):
					ax = self.lwpaxes[i_row, i_col]
					ax.entries = None
					ax.qns = None
					ax.corr_xs = ax.corr_ys = np.array([])
					threads.append(ax.update())

			thread.earlyreturn()

			for thread_ in threads:
				thread_.wait()
	
	@QThread.threaded_d
	@status_d
	@drawplot_decorator.d
	def set_data(self, thread=None):
		n_rows = config['plot_rows']
		n_cols = config['plot_cols']

		threads = []
		for asap_ax in self.lwpaxes.flatten():
			threads.append(asap_ax.update())

		for thread_ in threads:
				thread_.wait()
	
	def on_hover(self, event):
		x = event.xdata
		y = event.ydata

		if not all([x, y, event.inaxes]):
			text_cursor = ""
		else:
			if config['flag_showmainplotposition']:
				text_cursor = f"  ({x=:{config['flag_xformatfloat']}}, {y=:{config['flag_xformatfloat']}})  "
			else:
				text_cursor = ""
		mainwindow.statusbar.position_label.setText(text_cursor)
	
	def on_click(self, event):
		super().on_click(event)

		ax = event.inaxes
		x = event.xdata
		is_ctrl_pressed = (QApplication.keyboardModifiers() == Qt.KeyboardModifier.ControlModifier)

		if ax and x and is_ctrl_pressed:
			asap_ax = self.lwpaxes[self.__class__._active_ax_index]
			asap_ax.fit_data(x, x, onclick=True)

	def contextMenuCanvas(self, event):
		x, y = event.x(), event.y()
		geometry = self.plotcanvas.geometry()
		width, height = geometry.width(), geometry.height()
		x_rel, y_rel = x/width, 1 - y/height

		for lwpax in self.lwpaxes.flatten():
			xmin, ymin, width, height = lwpax.ax.get_position().bounds
			if xmin <= x_rel <= xmin + width and ymin <= y_rel <= ymin+height:
				break
		else: # Clicked outside of ax
			return

		menu = QMenu(self)
		set_active_action = menu.addAction('Make active Ax')
		fit_all_action = menu.addAction('Fit all')

		action = menu.exec(self.mapToGlobal(event.pos()))
		if action == set_active_action:
			mainwindow.lwpwidget._active_ax_index = (lwpax.row_i, lwpax.col_i)
		elif action == fit_all_action:
			i_col = lwpax.col_i
			AssignAllDialog.show_dialog(i_col)

class LASAPMainWindow(MainWindow):
	mainwidget_class = ASAPWidget
	menu_class = ASAPMenu

class LASAP(LLWP):
	mainwindow_class = LASAPMainWindow
	
	def debug_setup(self):
		pass



class ASAPDetailViewer(EQDockWidget):
	default_position = None
	available_in = ['LASAP',]
	
	drawplot = pyqtSignal()

	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		self.setWindowTitle("ASAP Detail Viewer")
	
		widget = QGroupBox()
		layout = QVBoxLayout()
		self.setWidget(widget)
		widget.setLayout(layout)
		
		tmp_layout = QHBoxLayout()
		tmp_layout.addWidget(QQ(QLabel, text='Width:'))
		tmp_layout.addWidget(QQ(QDoubleSpinBox, 'asap_detailviewerwidth'))
		tmp_layout.addStretch()

		layout.addLayout(tmp_layout)
		layout.addWidget(QQ(QCheckBox, 'asap_detailviewerfilter', text='Show only transitions used in cross-correlation'))

		self.fig = matplotlib.figure.Figure(dpi=config["plot_dpi"])
		self.plot_canvas = FigureCanvas(self.fig)
		layout.addWidget(self.plot_canvas)

		mainwindow.lwpwidget.peak_fitted.connect(self.update_view)
		self.drawplot.connect(self.draw_canvas)

	def update_view(self, asap_ax, offset):
		if not self.isVisible():
			return
			
		for ax in self.fig.get_axes():
			self.fig.delaxes(ax)
		
		# Prepare Entries
		entries = asap_ax.entries.copy()

		if config['asap_detailviewerfilter']:
			entries = entries[entries['use_for_cross_correlation']].copy()

		width = config['asap_detailviewerwidth'] or config['plot_width']
		entries['xmin'] = entries['x'] - width/2 + offset
		entries['xmax'] = entries['x'] + width/2 + offset
		entries['min_index'], entries['max_index'] = ExpFile.xs_to_indices(entries['xmin'], entries['xmax'])			
		
		gridspec_kw = {"hspace": 0, "wspace": 0}
		axes = self.fig.subplots(len(entries), gridspec_kw=gridspec_kw, squeeze=False)[::-1, 0]

		annotate_kwargs = {"x": 1, "y": 0.95, "horizontalalignment": "right", "verticalalignment": "top", "color": matplotlib.rcParams['text.color'], 'fontsize': 'small'}
		qns_labels = [[f'qn{ul}{i+1}' for i in range(config['series_qns'])] for ul in 'ul']

		for ax, (i, row) in zip(axes, entries.iterrows()):
			ax.xaxis.set_visible(False)
			ax.axvline(offset, color=config['color_fit'])

			min_index, max_index = row['min_index'], row['max_index']
			dataframe = ExpFile.df.iloc[min_index:max_index].copy()
			dataframe = dataframe[dataframe['visible']]
			
			if not len(dataframe):
				continue
			
			xs = dataframe['x'] - row['x']
			ys = dataframe['y']

			qnus_string = ','.join([f'{row[qn]}' for qn in qns_labels[0]])
			qnls_string = ','.join([f'{row[qn]}' for qn in qns_labels[1]])
			qnstring = f'{qnus_string} ← {qnls_string}'

			ax.plot(xs, ys, color=config['color_exp'])
			ax.text(**annotate_kwargs, s=qnstring, transform=ax.transAxes)
			ax.margins(0, 0.1)
		
		ax = axes[0]
		ax.xaxis.set_visible(True)
		ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(3))
		
		self.drawplot.emit()
	
	def draw_canvas(self):
		self.plot_canvas.draw_idle()

class ASAPSettingsWindow(ReferenceSeriesWindow):
	default_visible = True
	available_in = ['LASAP',]

	def __init__(self, *args, **kwargs):
		super(ReferenceSeriesWindow, self).__init__(*args, **kwargs)
		self.setWindowTitle("ASAP Settings")

		widget = QGroupBox()
		layout = QVBoxLayout(margin=True)
		self.setWidget(widget)
		widget.setLayout(layout)

		self.tab = QTabWidget()
		self.tab_order = None
		layout.addWidget(self.tab)

		self.tab.setTabsClosable(True)
		self.tab.setMovable(True)
		self.tab.setDocumentMode(True)

		self.tab.setTabBarAutoHide(True)
		self.tab.setCornerWidget(QQ(QToolButton, text="Dupl.", tooltip="Duplicate current tab", change=self.duplicate_tab), Qt.Corner.TopRightCorner)
		self.tab.tabCloseRequested.connect(self.close_tab)
		self.tab.tabBarDoubleClicked.connect(self.renameoradd_tab)
		self.tab.setCurrentIndex(config["series_currenttab"])

		self.tab_order = [self.tab.widget(i) for i in range(self.tab.count())]
		self.set_state(config["series_references"])
		if not self.tab.count():
			self.add_tab()

		mainwindow.lwpwidget.plotscreated.connect(self.update_number_of_references)
		self.tab.currentChanged.connect(self.check_order)
		
		self.toolbox = QToolBox()
		layout.addWidget(self.toolbox)

		settings_widget = QWidget()
		self.toolbox.addItem(settings_widget, 'Settings')
		tmp_layout = QGridLayout(margin=True)
		settings_widget.setLayout(tmp_layout)

		row_i = 0

		tmp_layout.addWidget(QQ(QLabel, text='Weighted Transitions: '), row_i, 0)
		tmp_layout.addWidget(QQ(QCheckBox, 'asap_weighted'), row_i, 1)

		row_i += 1

		tmp_layout.addWidget(QQ(QLabel, text='Interp. Resolution: '), row_i, 0)
		tmp_layout.addWidget(QQ(QDoubleSpinBox, 'asap_resolution', range=(0, None)), row_i, 1)

		row_i += 1

		tmp_layout.addWidget(QQ(QLabel, text='Only keep latest results: '), row_i, 0)
		tmp_layout.addWidget(QQ(QCheckBox, 'flag_keeponlylastassignment'), row_i, 1)
		
		tmp_layout.setRowStretch(row_i + 1, 1)


		egy_widget = QWidget()
		self.toolbox.addItem(egy_widget, 'Egy File')
		tmp_layout = QGridLayout(margin=True)
		egy_widget.setLayout(tmp_layout)

		row_i = 0

		tmp_layout.addWidget(QQ(QLabel, text='Energy File: '), row_i, 0)
		self.egy_file_button = QQ(QPushButton, text='Load File', change=lambda x: ASAPAx.load_egy_file())
		tmp_layout.addWidget(self.egy_file_button, row_i, 1)
		
		row_i += 1

		tmp_layout.addWidget(QQ(QLabel, text='Units Cat File:'), row_i, 0)
		tmp_layout.addWidget(QQ(QDoubleSpinBoxFullPrec, 'asap_catunitconversionfactor'), row_i, 1)

		
		tmp_layout.setRowStretch(row_i + 1, 1)


		filter_widget = QWidget()
		self.toolbox.addItem(filter_widget, 'Filter')
		tmp_layout = QVBoxLayout(margin=True)
		filter_widget.setLayout(tmp_layout)
		
		tmp_layout.addWidget(QQ(QLabel, text='Filter for transitions: '))
		tmp_layout.addWidget(QQ(QPlainTextEdit, 'asap_query'))

		tmp_layout.addStretch()

		layout.addWidget(QQ(QPushButton, text='Calculate Cross Correlation', change=lambda _: mainwindow.lwpwidget.calc_correlation_plots()))


	def add_tab(self, init_values={}, check_order=True):
		title = init_values.get("title", "Series")
		tmp = LevelSelector(self, init_values)
		tmp.values_changed.connect(self.changed)
		self.tab.addTab(tmp, title)
		if check_order:
			self.check_order()




##
## Startup
##

def start_llwp():
	global APP_TAG
	APP_TAG = 'LLWP'
	LLWP()

def start_lasap():
	File.special_file_handler = SpecialFilesHandlerASAP()
	AssignAllDialog.update_gui = AssignAllDialog.update_gui_asap

	global APP_TAG
	APP_TAG = 'LASAP'
	LASAP()

if __name__ == '__main__':
	# start_lasap()
	start_llwp()


##
## To Do
##

# - Auto detect format of .lin file (-> function in pyckett to check if current setting is reasonable for the file)
# - Which of these modules should be kept?
	# - EnergyLevelsTrendWindow
	# - SpectraResolverWindow
	# - CalibrateSpectrumWindow


##
## Some Tricks and Tips
##


## Open multiple files via glob string:

# import glob
# files = glob.glob(globstring)
# File.add_multiple_files_by_type(files)


## Hotkey to change series selector in specific way (increase Ka, decrease Kc)

# def tmp_function(change_value):
    # tab_widget = ReferenceSeriesWindow.instance.tab
    # refwidget = tab_widget.widget(config['series_currenttab'])
    # refwidget.setCurrentIndex(0)
    # seriesselector = refwidget.series_selector

    # noq = config['series_qns']
    # current_state = seriesselector.state
    # qnus = current_state['qnus'][:noq]
    # qnls = current_state['qnls'][:noq]
    
    # qnus[1] += change_value
    # qnus[2] -= change_value

    # qnls[1] += change_value
    # qnls[2] -= change_value
    
    # current_state['qnus'][:noq] = qnus
    # current_state['qnls'][:noq] = qnls
    # seriesselector.set_state()

# QShortcut('Ctrl+Y', mainwindow).activated.connect(lambda tmp_function=tmp_function: tmp_function(1))
# QShortcut('Ctrl+Shift+Y', mainwindow).activated.connect(lambda tmp_function=tmp_function: tmp_function(-1))

