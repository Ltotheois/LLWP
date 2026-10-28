#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Fit methods and lineshapes used by LLWP.

All methods inherit from :class:`FitMethod` and register themselves automatically
(see ``FitMethod.__init_subclass__``). The program lists the registered classes in
its menus, depending on which capabilities a class declares:

* ``can_fit_single``: fits a single line (main axes, assign all, ...)
* ``can_fit_multi``:  fits several lines at once (Blended Lines window)
* ``can_convolve``:   can create a kernel to convolute a stick spectrum

Hierarchy::

    FitMethod
    ├── Pgopher, Polynom, MultiPolynom    (determine a center position only)
    └── ParametricLineshape               (generic fit machinery for parametrized lines)
        ├── Gauss, Lorentz, Voigt
        └── FIDLineshape

A new lineshape only has to subclass :class:`ParametricLineshape`, declare its
parameters in ``params`` and implement ``profile``. The parameters are also used by
the Blended Lines window to create the columns of the results table.
"""

import numpy as np

from types import SimpleNamespace
from dataclasses import dataclass, field
from scipy import optimize, special

FWHM_TO_SIGMA = 1 / (2 * np.sqrt(2 * np.log(2)))


##
## Helper classes
##
@dataclass(frozen=True)
class Param:
    """Description of one parameter of a line.

    role:       "center", "amplitude" or "width" (used to create initial guesses and
                to find the position of a line) or None for any other parameter.
    shareable:  The value can be forced to be the same for all lines of a multi fit
                (e.g. the width).
    always_shared: The value is the same for all lines of a multi fit, independent of
                the user's choice (e.g. a property of the molecular beam).
    default, bounds: Initial guess and bounds for parameters without a role.
    first_line_value: If set, the parameter is fixed to this value for the first line
                and only fitted for the others (e.g. a phase relative to the first line).
    in_kernel:  The parameter without a role is also an input of the convolution kernel
                (parameters with the role "width" always are). Its `default` is the
                start value in the convolution window.
    """

    name: str
    label: str
    role: str = None
    shareable: bool = False
    always_shared: bool = False
    default: float = 0.0
    bounds: tuple = (-np.inf, np.inf)
    first_line_value: float = None
    in_kernel: bool = False


@dataclass
class Guess:
    """Information to create the initial guess and bounds of the parameters of one line."""

    index: int
    xmin: float
    xmax: float
    x: float
    y: float
    amp_min: float
    amp_max: float
    w0: float
    wmin: float
    wmax: float
    yscale: float = 1.0


@dataclass
class LastFit:
    label: str
    values: dict


@dataclass
class SingleFitResult:
    xmiddle: float
    xuncert: float
    fit_xs: np.ndarray
    fit_ys: np.ndarray
    popt: np.ndarray = None
    perr: np.ndarray = None
    values: dict = None  # Parameter values by name, only for parametric lineshapes


@dataclass
class MultiFitResult:
    """Result of fitting several lines at once.

    lines, errors: One dict per line, mapping the parameter names to their values.
    """

    method: "ParametricLineshape"
    lines: list
    errors: list
    baseline: list
    baseline_errors: list
    rms: float
    xcenter: float
    empty: bool = False
    context: SimpleNamespace = None
    raw_lines: list = field(default_factory=list, repr=False)
    raw_baseline: list = field(default_factory=list, repr=False)

    def evaluate(self, x):
        """Total fit curve including the baseline."""
        if self.empty:
            return np.zeros_like(x, dtype=float)
        y = self.method._evaluate(
            x, self.raw_lines, self.raw_baseline, self.context, self.xcenter
        )
        return y * self.context.yscale

    def evaluate_line(self, x, index):
        """Curve of a single line without the baseline."""
        if self.empty:
            return np.zeros_like(x, dtype=float)
        y = self.method.model(x, [self.raw_lines[index]], self.context)
        return y * self.context.yscale

    def evaluate_baseline(self, x):
        if not len(self.baseline):
            return np.zeros_like(x, dtype=float)
        return np.polyval(self.baseline, x - self.xcenter)

    def mean_values(self):
        """Mean value of each parameter over all lines."""
        if not self.lines or self.empty:
            return {}
        return {
            name: float(np.mean([line[name] for line in self.lines]))
            for name in self.lines[0]
        }

    def center(self, index):
        return self.lines[index][self.method.param_with_role("center").name]

    def amplitude(self, index):
        return self.lines[index][self.method.param_with_role("amplitude").name]


##
## Base class and registry
##
class FitMethod:
    """Base class of all fit methods. Subclasses are registered automatically."""

    registry = {}
    # Functions called with the class whenever a fit method is registered, e.g. to update menus
    registration_callbacks = []
    last_fit = None  # LastFit of the most recent fit the user did on an example line

    abstract = True  # Abstract classes are not registered
    name = None  # Name shown to the user
    derivatives = (0,)  # Supported derivative orders

    can_fit_single = False
    can_fit_multi = False
    can_convolve = False
    min_points = 2  # Minimum number of data points for a single fit

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        if cls.__dict__.get("abstract", False):
            return
        cls.abstract = False
        if cls.name is None:
            raise TypeError(f"The fit method {cls.__name__} needs a name.")
        FitMethod.registry[cls.name] = cls
        for callback in list(FitMethod.registration_callbacks):
            try:
                callback(cls)
            except Exception as E:
                print(f"Updating after registering {cls.name} failed: {E}")

    def __init__(self, derivative=0, **options):
        if derivative not in self.derivatives:
            raise NotImplementedError(
                f"The derivative {derivative} is not implemented for {self.name}."
            )
        self.derivative = derivative
        self.options = options

    def option(self, name, default=None):
        value = self.options.get(name)
        return default if value is None else value

    # Registry
    @staticmethod
    def make_label(name, derivative=0):
        if not derivative:
            return name
        ordinal = {1: "1st", 2: "2nd", 3: "3rd"}.get(derivative, f"{derivative}th")
        return f"{name} {ordinal} Derivative"

    @classmethod
    def labels(cls, capability=None):
        """Labels of all registered methods that support the capability, e.g. 'can_fit_multi'."""
        classes = [
            c for c in cls.registry.values() if not capability or getattr(c, capability)
        ]
        derivatives = sorted({d for c in classes for d in c.derivatives})
        return [
            cls.make_label(c.name, d)
            for d in derivatives
            for c in classes
            if d in c.derivatives
        ]

    @classmethod
    def resolve(cls, label):
        for method in cls.registry.values():
            for derivative in method.derivatives:
                if cls.make_label(method.name, derivative) == label:
                    return method, derivative
        raise NotImplementedError(f"The fit method '{label}' is not implemented.")

    @classmethod
    def create(cls, label, **options):
        method, derivative = cls.resolve(label)
        return method(derivative, **options)

    @classmethod
    def convolution_names(cls):
        return [c.name for c in cls.registry.values() if c.can_convolve]

    @property
    def label(self):
        return self.make_label(self.name, self.derivative)

    def remember(self, values):
        """Stores the parameter values of a fit to be reused, e.g. for the convolution."""
        if values:
            FitMethod.last_fit = LastFit(self.label, dict(values))

    # Interface
    def fit_single(self, xs, ys, peakdirection, fit_xs):
        raise NotImplementedError(f"{self.name} cannot fit single lines.")

    def fit_multi(self, xs, ys, peaks, xrange, **kwargs):
        raise NotImplementedError(f"{self.name} cannot fit multiple lines.")

    def kernel(self, xs, amplitude, **widths):
        raise NotImplementedError(f"{self.name} cannot create a convolution kernel.")


##
## Methods that only determine a center position
##
class CenterOnlyMethod(FitMethod):
    abstract = True
    can_fit_single = True

    @staticmethod
    def _result(xs, ys, peakdirection, fit_xs, ys_fit):
        xmiddle = (
            fit_xs[np.argmin(ys_fit)]
            if peakdirection < 0
            else fit_xs[np.argmax(ys_fit)]
        )
        return SingleFitResult(xmiddle, 0, fit_xs, ys_fit)


class Pgopher(CenterOnlyMethod):
    """Weighted center of mass of the upper half of the peak."""

    name = "Pgopher"
    min_points = 1

    def fit_single(self, xs, ys, peakdirection, fit_xs):
        ymin, ymax = np.min(ys), np.max(ys)

        if peakdirection < 0:
            cutoff = ymin + (ymax - ymin) / 2
            mask = ys <= cutoff
        else:
            cutoff = ymax - (ymax - ymin) / 2
            mask = ys >= cutoff

        fit_xs = xs[mask]
        fit_ys = ys[mask]

        # The weights are relative to the minimum, the shown curve stays on the data
        weights = fit_ys - ymin
        xmiddle = np.sum(fit_xs * weights) / np.sum(weights)
        return SingleFitResult(xmiddle, 0, fit_xs, fit_ys)


class Polynom(CenterOnlyMethod):
    """Polynomial of fixed rank, the maximum (minimum) is the center."""

    name = "Polynom"

    def rank(self, xs, ys):
        return self.option("polynomrank", 2)

    def fit_single(self, xs, ys, peakdirection, fit_xs):
        polynom = np.poly1d(np.polyfit(xs, ys, self.rank(xs, ys)))
        return self._result(xs, ys, peakdirection, fit_xs, polynom(fit_xs))


class MultiPolynom(Polynom):
    """Polynomial with the rank that yields the lowest residuals."""

    name = "MultiPolynom"

    def rank(self, xs, ys):
        maxrank = min(len(xs), self.option("polynommaxrank", 10))
        best_rms = np.inf
        best_rank = 0
        for rank in range(maxrank):
            polynom = np.poly1d(np.polyfit(xs, ys, rank))
            rms = np.mean((polynom(xs) - ys) ** 2)
            if rms < best_rms:
                best_rms = rms
                best_rank = rank
        return best_rank


##
## Parametrized lineshapes
##
class ParametricLineshape(FitMethod):
    """Base class for lines described by a number of parameters per line.

    Subclasses define ``params`` and ``profile`` (and optionally override ``model``,
    ``prepare`` and ``initial_guess`` for lineshapes that do not simply add up the
    profiles of the individual lines).
    """

    abstract = True
    can_fit_single = True
    can_fit_multi = True
    supports_baseline = True
    markers_at_amplitude = True  # Place the markers of the lines at their amplitude

    params = ()

    # Description of the lineshape
    def profile(self, x, *values):
        """Value of a single line, `values` are ordered like `params`."""
        raise NotImplementedError

    def model(self, x, lines, context):
        """Sum of the lines. `lines` is a list of value lists ordered like `params`."""
        ys = np.zeros_like(x, dtype=float)
        for line in lines:
            ys += self.profile(x, *line)
        return ys

    def prepare(self, xs, ys, xrange):
        """Creates the context passed to `model` and `initial_guess`."""
        xmin, xmax = xrange
        return SimpleNamespace(
            xmin=xmin, xmax=xmax, xcenter=(xmin + xmax) / 2, yscale=1.0
        )

    def initial_guess(self, guess):
        """List of (p0, lower bound, upper bound) for each parameter of a line."""
        result = []
        for param in self.params:
            if param.role == "center":
                result.append((guess.x, guess.xmin, guess.xmax))
            elif param.role == "amplitude":
                result.append((guess.y, guess.amp_min, guess.amp_max))
            elif param.role == "width":
                result.append((guess.w0, guess.wmin, guess.wmax))
            else:
                result.append((param.default, *param.bounds))
        return result

    # Convolution
    @classmethod
    def kernel_params(cls):
        return [p for p in cls.params if p.role == "width" or p.in_kernel]

    def kernel(self, xs, amplitude, **widths):
        values = []
        for param in self.params:
            if param.role == "center":
                values.append(0)
            elif param.role == "amplitude":
                values.append(amplitude)
            else:
                values.append(widths.get(param.name, param.default))
        return self.profile(xs, *values)

    # Parameter handling
    @classmethod
    def _names(cls):
        return [p.name for p in cls.params]

    @classmethod
    def param_with_role(cls, role):
        for param in cls.params:
            if param.role == role:
                return param
        raise KeyError(f"{cls.name} has no parameter with the role '{role}'.")

    @classmethod
    def _index_of_role(cls, role):
        return cls.params.index(cls.param_with_role(role))

    def _slots(self, nlines, shared):
        """Entries (line, parameter index) of the parameters that are varied during the fit.

        Parameters shared by all lines are listed once with line None after those of the lines.
        """
        slots = []
        for i in range(nlines):
            for k, param in enumerate(self.params):
                if self._is_shared(param, shared):
                    continue
                if i == 0 and param.first_line_value is not None:
                    continue
                slots.append((i, k))
        if nlines:
            slots.extend(
                (None, k)
                for k, p in enumerate(self.params)
                if self._is_shared(p, shared)
            )
        return slots

    @staticmethod
    def _is_shared(param, shared):
        return param.always_shared or (shared and param.shareable)

    def _unpack(self, vector, nlines, shared, nbaseline, fixed_value=None):
        """Splits the fit parameters into values of the lines and the baseline coefficients.

        fixed_value replaces the fixed values of parameters, e.g. to unpack errors."""
        lines = [[0.0] * len(self.params) for _ in range(nlines)]
        if nlines:
            for k, param in enumerate(self.params):
                if param.first_line_value is not None:
                    lines[0][k] = (
                        param.first_line_value if fixed_value is None else fixed_value
                    )

        slots = self._slots(nlines, shared)
        for value, (i, k) in zip(vector, slots):
            if i is None:
                for line in lines:
                    line[k] = value
            else:
                lines[i][k] = value
        baseline = list(vector[len(slots) : len(slots) + nbaseline])
        return lines, baseline

    def _rescale(self, vector, nlines, shared, nbaseline, scale):
        """Converts amplitudes and baseline between the fit and the data scale."""
        vector = np.array(vector, dtype=float)
        slots = self._slots(nlines, shared)
        for j, (_, k) in enumerate(slots):
            if self.params[k].role == "amplitude":
                vector[j] *= scale
        vector[len(slots) : len(slots) + nbaseline] *= scale
        return vector

    def _evaluate(self, x, lines, baseline, context, xcenter):
        ys = self.model(x, lines, context)
        if len(baseline):
            ys = ys + np.polyval(baseline, x - xcenter)
        return ys

    def _curve_fit(
        self, xs, ys, context, nlines, shared, nbaseline, xcenter, p0, bounds
    ):
        def function(x, *vector):
            lines, baseline = self._unpack(vector, nlines, shared, nbaseline)
            return self._evaluate(x, lines, baseline, context, xcenter)

        popt, pcov = optimize.curve_fit(function, xs, ys, p0=p0, bounds=bounds)
        perr = np.sqrt(np.diag(pcov))
        return function, popt, perr

    def _slot_guesses(self, guesses, slots):
        p0, lower, upper = [], [], []
        for i, k in slots:
            value, lo, hi = guesses[i or 0][k]
            p0.append(value)
            lower.append(lo)
            upper.append(hi)
        return p0, lower, upper

    # Single line fit
    def fit_single(self, xs, ys, peakdirection, fit_xs):
        xmin, xmax = xs.min(), xs.max()
        x0 = (xmin + xmax) / 2

        xs_weight_factor = self.option("xs_weight_factor", 4)
        if xs_weight_factor:
            ys_weighted = ys * np.exp(
                -np.abs(np.abs(xs - x0) / (xmax - xmin)) * xs_weight_factor
            )
            x0 = (
                xs[np.argmax(ys_weighted)]
                if peakdirection >= 0
                else xs[np.argmin(ys_weighted)]
            )

        ymin, ymax, ymean, yptp = ys.min(), ys.max(), ys.mean(), np.ptp(ys)
        y0 = 0

        w0 = self.option("w0", (xmax - xmin) / 10)
        wmin = self.option("wmin", 0)
        wmax = self.option("wmax", (xmax - xmin))

        amp_min, amp_max = -3 * yptp, 3 * yptp
        if peakdirection < 0:
            amp_max = 0
            y0 = -yptp
        if peakdirection > 0:
            amp_min = 0
            y0 = yptp

        context = self.prepare(xs, ys, (xmin, xmax))
        scale = context.yscale
        guess = Guess(0, xmin, xmax, x0, y0, amp_min, amp_max, w0, wmin, wmax, scale)
        guesses = [self.initial_guess(guess)]

        slots = self._slots(1, False)
        p0, lower, upper = self._slot_guesses(guesses, slots)

        offset = bool(self.option("offset", False))
        if offset:
            p0.append(ymean / scale)
            lower.append(ymin / scale)
            upper.append(ymax / scale)

        function, popt, perr = self._curve_fit(
            xs, ys / scale, context, 1, False, int(offset), 0, p0, [lower, upper]
        )
        fit_ys = function(fit_xs, *popt) * scale

        popt = self._rescale(popt, 1, False, int(offset), scale)
        perr = self._rescale(perr, 1, False, int(offset), scale)
        k = self._index_of_role("center")
        xmiddle = self._unpack(popt, 1, False, 0)[0][0][k]
        xuncert = self._unpack(perr, 1, False, 0, fixed_value=0)[0][0][k]

        values = dict(zip(self._names(), self._unpack(popt, 1, False, 0)[0][0]))
        return SingleFitResult(xmiddle, xuncert, fit_xs, fit_ys, popt, perr, values)

    # Multi line fit
    def fit_multi(
        self,
        xs,
        ys,
        peaks,
        xrange,
        baseline_rank=-1,
        shared=False,
        wmax=None,
        peakdirection=0,
    ):
        """Fits all peaks at once.

        peaks:         list of (x, y, x_rel) with the positions chosen by the user and
                       their position relative to the center of the range.
        baseline_rank: Rank of the baseline polynomial, -1 for no baseline.
        shared:        Use the same value for all parameters that are `shareable`.
        wmax:          Upper bound of the width parameters.
        """
        xmin, xmax = xrange
        xwidth = xmax - xmin
        xcenter = (xmin + xmax) / 2
        nlines = len(peaks)
        nbaseline = baseline_rank + 1 if self.supports_baseline else 0
        shared = bool(shared) and nlines > 0

        context = self.prepare(xs, ys, xrange)
        scale = context.yscale

        if not (len(ys) and (nbaseline + nlines)):
            return MultiFitResult(
                self,
                [{p.name: 0 for p in self.params} for _ in range(nlines)],
                [{p.name: 0 for p in self.params} for _ in range(nlines)],
                [0] * nbaseline,
                [0] * nbaseline,
                0,
                xcenter,
                empty=True,
                context=context,
            )

        yptp = 4 * (np.amax(ys) - np.amin(ys))
        w0 = xwidth
        wmax = wmax or w0
        amp_min, amp_max = -3 * yptp, 3 * yptp
        if peakdirection < 0:
            amp_max = 0
        if peakdirection > 0:
            amp_min = 0

        guesses = []
        for i, (x, y, x_rel) in enumerate(peaks):
            if not xmin < x < xmax:
                x = xcenter + x_rel
                if not xmin < x < xmax:
                    x = xcenter
                y = yptp * np.sign(peakdirection)
            elif not amp_min < y < amp_max:
                y = yptp * np.sign(peakdirection)

            guess = Guess(
                i, xmin, xmax, x, y, amp_min, amp_max, wmax / 10, 0, wmax, scale
            )
            guesses.append(self.initial_guess(guess))

        slots = self._slots(nlines, shared)
        p0, lower, upper = self._slot_guesses(guesses, slots)
        p0.extend([0] * nbaseline)
        lower.extend([-np.inf] * nbaseline)
        upper.extend([+np.inf] * nbaseline)

        function, popt, perr = self._curve_fit(
            xs,
            ys / scale,
            context,
            nlines,
            shared,
            nbaseline,
            xcenter,
            p0,
            [lower, upper],
        )
        raw_lines, raw_baseline = self._unpack(popt, nlines, shared, nbaseline)
        res_ys = function(xs, *popt) * scale
        rms = np.sqrt(np.sum((ys - res_ys) ** 2) / len(ys))

        popt = self._rescale(popt, nlines, shared, nbaseline, scale)
        perr = self._rescale(perr, nlines, shared, nbaseline, scale)
        lines, baseline = self._unpack(popt, nlines, shared, nbaseline)
        errors, baseline_errors = self._unpack(
            perr, nlines, shared, nbaseline, fixed_value=0
        )

        names = self._names()
        return MultiFitResult(
            self,
            [dict(zip(names, line)) for line in lines],
            [dict(zip(names, line)) for line in errors],
            baseline,
            baseline_errors,
            rms,
            xcenter,
            context=context,
            raw_lines=raw_lines,
            raw_baseline=raw_baseline,
        )


##
## Implementations of lineshapes
##
def gaussian(derivative, x, x0, amp, fwhm):
    sigma = fwhm * FWHM_TO_SIGMA
    if sigma == 0:
        return np.where(x == x0, np.inf, 0.0)

    if derivative == 0:
        ys = amp * np.exp(-((x - x0) ** 2) / (2 * sigma**2))
    elif derivative == 1:
        ys = (
            -amp
            / (sigma * np.exp(-0.5))
            * (x - x0)
            * np.exp(-((x - x0) ** 2) / (2 * sigma**2))
        )
    elif derivative == 2:
        ys = (
            amp
            * (1 - ((x - x0) / sigma) ** 2)
            * np.exp(-((x - x0) ** 2) / (2 * sigma**2))
        )
    else:
        raise NotImplementedError(
            "Only the zeroth, first, and second derivatives of a Gaussian are implemented."
        )
    return ys


def lorentzian(derivative, x, x0, amp, fwhm):
    gamma = fwhm / 2
    if gamma == 0:
        return np.where(x == x0, np.inf, 0.0)

    if derivative == 0:
        ys = amp * gamma**2 / (gamma**2 + (x - x0) ** 2)
    elif derivative == 1:
        ys = (
            (-amp * gamma**3 * 16 / 9 * np.sqrt(3))
            * (x - x0)
            / ((x - x0) ** 2 + gamma**2) ** 2
        )
    elif derivative == 2:
        ys = (
            (amp * gamma**4)
            * (gamma**2 - 3 * (x - x0) ** 2)
            / ((x - x0) ** 2 + gamma**2) ** 3
        )
    else:
        raise NotImplementedError(
            "Only the zeroth, first, and second derivatives of a Lorentzian are implemented."
        )
    return ys


def voigt(derivative, x, x0, amp, fwhm_gauss, fwhm_lorentz):
    sigma = fwhm_gauss * FWHM_TO_SIGMA
    gamma = fwhm_lorentz / 2
    if gamma == sigma == 0:
        return np.where(x == x0, np.inf, 0.0)

    z = (x - x0 + 1j * gamma) / (sigma * np.sqrt(2))
    wz = special.wofz(z)
    w0 = special.wofz((1j * gamma) / (sigma * np.sqrt(2)))

    if derivative == 0:

        def tmp(x, x0, wz, sigma, gamma):
            return np.real(wz) / (sigma * np.sqrt(2 * np.pi))
    elif derivative == 1:

        def tmp(x, x0, wz, sigma, gamma):
            return (
                1
                / (sigma**3 * np.sqrt(2 * np.pi))
                * (gamma * np.imag(wz) - (x - x0) * np.real(wz))
            )
    elif derivative == 2:

        def tmp(x, x0, wz, sigma, gamma):
            return (
                1
                / (sigma**5 * np.sqrt(2 * np.pi))
                * (
                    gamma * (2 * (x - x0) * np.imag(wz) - sigma * np.sqrt(2 / np.pi))
                    + (gamma**2 + sigma**2 - (x - x0) ** 2) * np.real(wz)
                )
            )
    else:
        raise NotImplementedError(
            "Only the zeroth, first, and second derivatives of a Voigt profile are implemented."
        )

    ys = tmp(x, x0, wz, sigma, gamma)
    if derivative == 1:
        xspan = max(sigma, gamma) * 2
        tmp_xs = np.linspace(-xspan, +xspan, 1000)
        ymax = np.max(tmp(tmp_xs, 0, w0, sigma, gamma))
    else:
        ymax = tmp(0, 0, w0, sigma, gamma)
    if ymax == 0:
        ymax = 1
    ys *= amp / ymax

    return ys


CENTER = Param("x0", "Frequency", role="center")
AMPLITUDE = Param("amp", "Amplitude", role="amplitude")
FWHM_GAUSS = Param("fwhm_gauss", "FWHM Gauss", role="width", shareable=True)
FWHM_LORENTZ = Param("fwhm_lorentz", "FWHM Lorentz", role="width", shareable=True)


class Gauss(ParametricLineshape):
    name = "Gauss"
    derivatives = (0, 1, 2)
    can_convolve = True
    params = (CENTER, AMPLITUDE, FWHM_GAUSS)

    def profile(self, x, x0, amp, fwhm):
        return gaussian(self.derivative, x, x0, amp, fwhm)


class Lorentz(ParametricLineshape):
    name = "Lorentz"
    derivatives = (0, 1, 2)
    can_convolve = True
    params = (CENTER, AMPLITUDE, FWHM_LORENTZ)

    def profile(self, x, x0, amp, fwhm):
        return lorentzian(self.derivative, x, x0, amp, fwhm)


class Voigt(ParametricLineshape):
    name = "Voigt"
    derivatives = (0, 1, 2)
    can_convolve = True
    params = (CENTER, AMPLITUDE, FWHM_GAUSS, FWHM_LORENTZ)

    def profile(self, x, x0, amp, fwhm_gauss, fwhm_lorentz):
        return voigt(self.derivative, x, x0, amp, fwhm_gauss, fwhm_lorentz)


class FIDLineshape(ParametricLineshape):
    """Fits the power spectrum of a sum of exponentially decaying oscillations.

    The lines are added in the time domain, therefore their phases matter. The phase of
    the first line is fixed to zero. To process the FID before the Fourier transform set
    FIDLineshape.postprocessing_function to a function (ts, ys_fid) -> (ts, ys_fid).
    """

    name = "FID Fit"
    supports_baseline = False
    markers_at_amplitude = False
    postprocessing_function = None

    params = (
        CENTER,
        AMPLITUDE,
        Param("fwhm", "FWHM", role="width", shareable=True),
        Param("phase", "Phase", default=0, bounds=(-np.pi, np.pi), first_line_value=0),
    )

    def prepare(self, xs, ys, xrange):
        context = super().prepare(xs, ys, xrange)
        xwidth = context.xmax - context.xmin
        resolution = self.option("fid_resolution", 0.05)

        yscale = np.max(ys) if len(ys) else 1
        context.yscale = yscale if yscale > 0 else 1
        context.dt = 1 / xwidth
        context.n_time = max(int(1 / context.dt / resolution), 2)
        context.ts = np.arange(context.n_time) * context.dt
        return context

    def initial_guess(self, guess):
        xwidth = guess.xmax - guess.xmin
        ymax = 10
        y = guess.y / guess.yscale
        if not 0 < y < ymax:
            y = ymax / 2
        w0 = min(guess.wmax / 2, xwidth / 2)
        phase = self.params[3]
        return [
            (guess.x, guess.xmin, guess.xmax),
            (y, 0, ymax),
            (w0, 0, guess.wmax),
            (phase.default, *phase.bounds),
        ]

    def line_signal(self, ts, context, x0, y0, fwhm0, phi0):
        """Time signal of a single line."""
        return (
            y0
            * np.pi
            * fwhm0
            * np.exp(-ts * np.pi * fwhm0)
            * np.exp(1j * (2 * np.pi * (x0 - context.xcenter) * ts + phi0))
        )

    def model(self, x, lines, context):
        ts = context.ts
        ys_fid = np.zeros_like(ts, dtype=complex)
        for line in lines:
            ys_fid += self.line_signal(ts, context, *line)

        # Accessed through the class so that plain functions are not bound to the instance
        postprocessing = type(self).postprocessing_function
        if postprocessing is not None:
            ts, ys_fid = postprocessing(ts, ys_fid)

        spec = np.fft.fftshift(np.fft.fft(ys_fid, context.n_time)) * context.dt
        freq = np.fft.fftshift(np.fft.fftfreq(context.n_time, context.dt))
        freq += context.xcenter
        return np.interp(x, freq, np.abs(spec))
