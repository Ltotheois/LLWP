"""Helpers for the pyqtgraph plots of LLWP.

The data for the plots is calculated in worker threads, while the graphics items of
pyqtgraph must only be touched in the GUI thread. The plots therefore describe what to
draw in plain data (the *spec*) and pass it to :func:`run_in_gui` together with a
function that applies it to the items.

A spec of a spectrum has the following keys:

* ``exp`` and ``cat``: lists of ``(xs, ys, connect, color)``, one for each curve
  (``connect`` is ``"all"`` for a line or ``"pairs"`` for sticks, see :func:`stick_xy`)
* ``lin``: ``(xs, ys, colors)`` of the markers of the assignments
* ``yrange``: ``(ymin, ymax)``
"""

import os
import traceback as tb
from xml.sax.saxutils import escape as xml_escape

import numpy as np
import pandas as pd
import pyqtgraph as pg
from pyqtgraph.graphicsItems.ScatterPlotItem import Symbols

from PyQt6.QtCore import QObject, QPointF, QRectF, Qt, pyqtSignal
from PyQt6.QtGui import (
    QColor,
    QFont,
    QFontInfo,
    QFontMetricsF,
    QGuiApplication,
    QPainter,
    QPainterPath,
    QPalette,
    QPen,
    QPicture,
)
from PyQt6.QtWidgets import QApplication, QFileDialog, QGraphicsRectItem, QMenu


##
## Threading
##
class GuiInvoker(QObject):
    """Runs functions in the thread in which this object was created (the GUI thread)."""

    requested = pyqtSignal(object)

    def __init__(self):
        super().__init__()
        self.requested.connect(self._run)

    def _run(self, func):
        try:
            func()
        except Exception:
            tb.print_exc()


_gui_invoker = None


def run_in_gui(func):
    """Calls `func` in the GUI thread (directly, if called from there).

    The first call has to be made from the GUI thread."""
    global _gui_invoker
    if _gui_invoker is None:
        _gui_invoker = GuiInvoker()
    _gui_invoker.requested.emit(func)


##
## Theme
##
class ThemeNotifier(QObject):
    """Emits `changed` when the colors of the plots have to be updated and
    `style_changed` when other things of their look changed (the rounded corners)."""

    changed = pyqtSignal()
    style_changed = pyqtSignal()


theme = ThemeNotifier()


def is_dark_theme():
    return QGuiApplication.instance().styleHints().colorScheme() == Qt.ColorScheme.Dark


def plot_foreground():
    return "#ffffff" if is_dark_theme() else "#000000"


def plot_background():
    return "#000000" if is_dark_theme() else "#ffffff"


def update_pyqtgraph_options():
    """Sets the options of pyqtgraph that all plots share, they follow the theme."""
    pg.setConfigOptions(
        foreground=plot_foreground(), background=plot_background(), antialias=True
    )


theme.changed.connect(update_pyqtgraph_options)


def blend(first, second, fraction):
    """The color `fraction` of the way from `first` to `second`."""
    first, second = QColor(first), QColor(second)
    return QColor(
        *(
            round(a + fraction * (b - a))
            for a, b in zip(first.getRgb()[:3], second.getRgb()[:3], strict=True)
        )
    )


def plot_frame_color():
    """The color of the frame of the plots and of their ticks: softer than the text."""
    return blend(plot_background(), plot_foreground(), 0.45)


def plot_inner_frame_color():
    """The gray of the frame of the group box around the plots."""
    return QApplication.palette().color(QPalette.ColorRole.Window).darker(111)


FONT_SIZE = 10  # in points


def plot_font(size=FONT_SIZE):
    font = QFont()
    font.setPointSizeF(float(size))
    return font


##
## Colors, pens and markers
##
# Colors that matplotlib knows as single letters, which are still allowed in the options
_SINGLE_LETTER_COLORS = {
    "b": "#0000ff",
    "g": "#008000",
    "r": "#ff0000",
    "c": "#00bfbf",
    "m": "#bf00bf",
    "y": "#bfbf00",
    "k": "#000000",
    "w": "#ffffff",
}


def parse_color(color):
    """QColor of a color name or of a hex color (#rgb, #rrggbb or #rrggbbaa).

    Contrary to Qt, the opacity of eight digit hex colors comes last."""
    color = _SINGLE_LETTER_COLORS.get(color, color)
    if isinstance(color, str) and color.startswith("#") and len(color) == 9:
        r, g, b, a = (int(color[i : i + 2], 16) for i in (1, 3, 5, 7))
        return QColor(r, g, b, a)

    qcolor = QColor(color)
    if not qcolor.isValid():
        raise ValueError(f"Invalid color: '{color}'")
    return qcolor


def to_hex(color):
    """The color as #rrggbb without opacity."""
    return parse_color(color).name()


def to_qcolor(color, alpha=1.0):
    qcolor = parse_color(color)
    qcolor.setAlphaF(qcolor.alphaF() * alpha)
    return qcolor


def make_pen(color, alpha=1.0, width=1, style=None):
    pen = pg.mkPen(to_qcolor(color, alpha), width=width)
    if style is not None:
        pen.setStyle(style)
    return pen


def is_symbol(symbol):
    """Is `symbol` the name of a symbol that pyqtgraph can draw."""
    return symbol in Symbols


def stick_xy(xs, ys, base=0):
    """Coordinates for vertical lines from (x, base) to (x, y) drawn with connect='pairs'."""
    x = np.repeat(np.asarray(xs, dtype=float), 2)
    y = np.empty(len(x))
    y[0::2] = base
    y[1::2] = ys
    return x, y


def split_by_color(colors):
    """Yields each distinct color and the boolean mask of its entries."""
    codes, uniques = pd.factorize(np.asarray(colors, dtype=object))
    for code, color in enumerate(uniques):
        yield color, codes == code


##
## Rounded corners
##
STYLE = {"radius": 4}  # in pixels, `plot_radius` of the options


def set_plot_radius(radius):
    """Sets the radius of the rounded corners of the plots in pixels (0 for square corners)."""
    STYLE["radius"] = max(int(radius), 0)
    theme.style_changed.emit()


class RoundedCornersMixin:
    """Paints rounded corners on a plot widget, in the color of the surroundings.

    The plots have a margin of at least `max_radius` pixels to the edge of the widget, so
    that the corners stay clear of the frame and the labels."""

    max_radius = 6
    is_square = False

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        theme.changed.connect(self.viewport().update)
        theme.style_changed.connect(self.viewport().update)

    def paintEvent(self, event):
        super().paintEvent(event)
        radius = 0 if self.is_square else min(STYLE["radius"], self.max_radius)
        view = self.viewport()
        if radius <= 0 or view.width() < 4 * radius or view.height() < 4 * radius:
            return

        outside, inside = QPainterPath(), QPainterPath()
        outside.addRect(QRectF(view.rect()))
        inside.addRoundedRect(QRectF(view.rect()), radius, radius)
        painter = QPainter(view)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.fillPath(
            outside.subtracted(inside),
            QApplication.palette().color(QPalette.ColorRole.Window),
        )
        painter.end()

    def grab_clean(self):
        """The plot as an image, with square corners."""
        self.is_square = True
        try:
            return self.grab()
        finally:
            self.is_square = False


class RoundedPlotWidget(RoundedCornersMixin, pg.PlotWidget):
    """A pg.PlotWidget with rounded corners, the plot has a margin of `margin` pixels."""

    margin = 10
    max_radius = 10

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.getPlotItem().layout.setContentsMargins(*[self.margin] * 4)


class RoundedGraphicsLayoutWidget(RoundedCornersMixin, pg.GraphicsLayoutWidget):
    """A pg.GraphicsLayoutWidget with rounded corners and margins for the plots."""

    max_radius = 6

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Wide enough for the labels of the ticks at the edges of the plots
        self.ci.layout.setContentsMargins(25, 6, 25, 6)
        self.ci.layout.setSpacing(0)


##
## Saving
##
##
## SVG export
##
# Qt can write an SVG of the widget, but pyqtgraph draws the lines, ticks and markers in
# a way that does not survive it (misplaced ticks, wrong lines). pyqtgraph's own exporter
# fails for the items used here. Therefore the SVG is written from the data of the items.
class SvgWriter:
    """Collects the elements of an SVG, with the clipping of the plots."""

    def __init__(self, width, height):
        self.width, self.height = width, height
        self.elements = []

    @staticmethod
    def color(qcolor):
        return f"rgb({qcolor.red()},{qcolor.green()},{qcolor.blue()})"

    @classmethod
    def stroke(cls, pen, unfilled=True):
        """Attributes of a pen, or None if it does not draw."""
        if pen is None or pen.style() == Qt.PenStyle.NoPen:
            return None
        color = pen.color()
        attributes = (
            f'stroke="{cls.color(color)}" stroke-opacity="{color.alphaF():.3f}" '
            f'stroke-width="{max(pen.widthF(), 1):g}"'
        )
        if unfilled:
            attributes += ' fill="none"'
        if pen.style() == Qt.PenStyle.DashLine:
            attributes += ' stroke-dasharray="4,3"'
        return attributes

    @classmethod
    def fill(cls, brush):
        if brush is None or brush.style() == Qt.BrushStyle.NoBrush:
            return 'fill="none"'
        color = brush.color()
        return f'fill="{cls.color(color)}" fill-opacity="{color.alphaF():.3f}"'

    def add(self, element):
        self.elements.append(element)

    def line(self, pen, p1, p2):
        stroke = self.stroke(pen)
        if stroke:
            self.add(
                f'<line x1="{p1.x():.2f}" y1="{p1.y():.2f}" '
                f'x2="{p2.x():.2f}" y2="{p2.y():.2f}" {stroke}/>'
            )

    def text(self, text, x, y, font, color, anchor="middle", vertical="middle"):
        """Text at x (`anchor` is start, middle or end) and y (`vertical` is which part
        of the text is at y: top, middle or bottom).

        The baseline of the text is calculated as not all viewers support the
        alignment of the baseline."""
        metrics = QFontMetricsF(font)
        baseline = {
            "top": y + metrics.ascent(),
            "middle": y + (metrics.ascent() - metrics.descent()) / 2,
            "bottom": y - metrics.descent(),
        }[vertical]
        size = QFontInfo(font).pixelSize()
        self.add(
            f'<text x="{x:.2f}" y="{baseline:.2f}" font-family="{font.family()}" '
            f'font-size="{size}" fill="{self.color(color)}" '
            f'text-anchor="{anchor}">{xml_escape(text)}</text>'
        )

    def tostring(self, background):
        return (
            '<?xml version="1.0" encoding="UTF-8"?>\n'
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{self.width}" '
            f'height="{self.height}" viewBox="0 0 {self.width} {self.height}">\n'
            f'<rect width="{self.width}" height="{self.height}" '
            f"{self.fill(background)}/>\n" + "\n".join(self.elements) + "\n</svg>\n"
        )


def _map(transform, xs, ys):
    """Maps the coordinates with a QTransform."""
    xs, ys = np.asarray(xs, dtype=float), np.asarray(ys, dtype=float)
    return (
        transform.m11() * xs + transform.m21() * ys + transform.dx(),
        transform.m12() * xs + transform.m22() * ys + transform.dy(),
    )


def _clip_segments(xs, ys, rect):
    """Clips the segments between the points to the rectangle (Liang-Barsky).

    Returns the start and end of each segment, which is None where it is outside, and
    the parts t0 and t1 of the segment that are inside."""
    xmin, xmax, ymin, ymax = rect.left(), rect.right(), rect.top(), rect.bottom()
    x0, y0, x1, y1 = xs[:-1], ys[:-1], xs[1:], ys[1:]
    dx, dy = x1 - x0, y1 - y0

    t0, t1 = np.zeros(len(dx)), np.ones(len(dx))
    valid = np.isfinite(x0) & np.isfinite(y0) & np.isfinite(x1) & np.isfinite(y1)
    for p, q in (
        (-dx, x0 - xmin),
        (dx, xmax - x0),
        (-dy, y0 - ymin),
        (dy, ymax - y0),
    ):
        parallel = p == 0
        valid &= ~(parallel & (q < 0))
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(parallel, 0.0, q / p)
        t0 = np.where(~parallel & (p < 0), np.maximum(t0, ratio), t0)
        t1 = np.where(~parallel & (p > 0), np.minimum(t1, ratio), t1)
    valid &= t0 <= t1
    return valid, t0, t1, (x0, y0, dx, dy)


def _clipped_polylines(xs, ys, rect):
    """The parts of the line through the points that are inside of the rectangle."""
    valid, t0, t1, (x0, y0, dx, dy) = _clip_segments(xs, ys, rect)
    indices = np.flatnonzero(valid)
    if not len(indices):
        return []

    # Segments that continue each other are one polyline
    continues = (
        (np.diff(indices) == 1) & (t1[indices[:-1]] == 1) & (t0[indices[1:]] == 0)
    )
    starts = np.flatnonzero(np.concatenate(([True], ~continues)))
    ends = np.concatenate((starts[1:], [len(indices)]))

    polylines = []
    for first, last in zip(starts, ends):
        group = indices[first:last]
        px = np.append(
            x0[group] + t0[group] * dx[group],
            x0[group[-1]] + t1[group[-1]] * dx[group[-1]],
        )
        py = np.append(
            y0[group] + t0[group] * dy[group],
            y0[group[-1]] + t1[group[-1]] * dy[group[-1]],
        )
        polylines.append((px, py))
    return polylines


def _svg_curve(svg, curve, transform, rect):
    xs, ys = curve.getData()
    if xs is None or not len(xs):
        return
    stroke = svg.stroke(curve.opts["pen"])
    if not stroke:
        return

    xs, ys = _map(transform, xs, ys)
    if curve.opts["connect"] == "pairs":
        n = len(xs) // 2
        # Each pair is a segment, they do not continue each other
        pair_xs = np.column_stack((xs[0 : 2 * n : 2], xs[1 : 2 * n : 2]))
        pair_ys = np.column_stack((ys[0 : 2 * n : 2], ys[1 : 2 * n : 2]))
        parts = []
        for px, py in zip(pair_xs, pair_ys):
            parts.extend(_clipped_polylines(px, py, rect))
    else:
        parts = _clipped_polylines(xs, ys, rect)

    path = " ".join(
        "M" + " L".join(f"{x:.2f},{y:.2f}" for x, y in zip(px, py)) for px, py in parts
    )
    if path:
        svg.add(f'<path d="{path}" {stroke}/>')


def _svg_marker(svg, symbol, x, y, size, pen, brush):
    """One marker of a scatter plot."""
    stroke = svg.stroke(pen, unfilled=False)
    style = f"{svg.fill(brush)} {stroke or ''}"
    r = size / 2
    if symbol in ("+", "x"):
        line_pen = pen if stroke else QPen(brush.color())
        if symbol == "+":
            ends = [((x - r, y), (x + r, y)), ((x, y - r), (x, y + r))]
        else:
            ends = [((x - r, y - r), (x + r, y + r)), ((x - r, y + r), (x + r, y - r))]
        for p1, p2 in ends:
            svg.line(line_pen, QPointF(*p1), QPointF(*p2))
    elif symbol == "s":
        svg.add(
            f'<rect x="{x - r:.2f}" y="{y - r:.2f}" width="{size:.2f}" '
            f'height="{size:.2f}" {style}/>'
        )
    elif symbol in ("t", "t1", "t2", "t3", "d", "star", "p", "h"):
        if symbol == "d":
            corners = [(0, -1), (1, 0), (0, 1), (-1, 0)]
        elif symbol == "star":
            corners = [
                (
                    np.sin(k * np.pi / 5) * (1 if k % 2 == 0 else 0.4),
                    -np.cos(k * np.pi / 5) * (1 if k % 2 == 0 else 0.4),
                )
                for k in range(10)
            ]
        elif symbol in ("p", "h"):
            n = 5 if symbol == "p" else 6
            corners = [
                (np.sin(2 * np.pi * k / n), -np.cos(2 * np.pi * k / n))
                for k in range(n)
            ]
        else:
            angle = {"t": np.pi, "t1": 0, "t2": np.pi / 2, "t3": -np.pi / 2}[symbol]
            corners = [
                (np.sin(angle + 2 * np.pi * k / 3), -np.cos(angle + 2 * np.pi * k / 3))
                for k in range(3)
            ]
        points = " ".join(f"{x + r * cx:.2f},{y + r * cy:.2f}" for cx, cy in corners)
        svg.add(f'<polygon points="{points}" {style}/>')
    else:
        svg.add(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{r:.2f}" {style}/>')


def _svg_scatter(svg, scatter, transform, rect):
    data = scatter.data
    if not len(data):
        return
    xs, ys = _map(transform, data["x"], data["y"])
    inside = (
        (xs >= rect.left())
        & (xs <= rect.right())
        & (ys >= rect.top())
        & (ys <= rect.bottom())
    )
    for x, y, size, symbol, pen, brush in zip(
        xs[inside],
        ys[inside],
        data["size"][inside],
        data["symbol"][inside],
        data["pen"][inside],
        data["brush"][inside],
    ):
        # The default of the item is used where a point has no own style (size -1, ...)
        size = scatter.opts["size"] if size < 0 else size
        symbol = scatter.opts["symbol"] if symbol is None else symbol
        pen = scatter.opts["pen"] if pen is None else pen
        brush = scatter.opts["brush"] if brush is None else brush
        if isinstance(symbol, str):
            _svg_marker(svg, symbol, x, y, size, pen, brush)


def _svg_text(svg, item, transform):
    text = item.textItem.toPlainText()
    if not text:
        return
    # The text keeps its size, so only its position is transformed
    position = transform.map(QPointF(0, 0))
    font = item.textItem.font()
    rect = item.textItem.boundingRect()  # with the margin around the text
    margin = item.textItem.document().documentMargin()
    left = position.x() - item.anchor[0] * rect.width() + margin
    top = position.y() - item.anchor[1] * rect.height() + margin
    line_height = QFontMetricsF(font).lineSpacing()
    for i, line in enumerate(text.split("\n")):
        svg.text(
            line,
            left,
            top + i * line_height,
            font,
            item.textItem.defaultTextColor(),
            anchor="start",
            vertical="top",
        )


def _svg_axis(svg, axis, transform):
    picture = QPicture()
    painter = QPainter(picture)
    try:
        specs = axis.generateDrawSpecs(painter)
    finally:
        painter.end()
    if specs is None:
        return
    axis_spec, tick_specs, text_specs = specs

    to_scene = axis.sceneTransform() * transform

    def mapped(point):
        x, y = _map(to_scene, [point.x()], [point.y()])
        return QPointF(x[0], y[0])

    if axis_spec is not None:
        pen, p1, p2 = axis_spec
        svg.line(pen, mapped(p1), mapped(p2))
    for pen, p1, p2 in tick_specs:
        svg.line(pen, mapped(p1), mapped(p2))

    font = axis.style.get("tickFont") or axis.font()
    color = axis.textPen().color()
    for rect, flags, string in text_specs:
        horizontal = flags & Qt.AlignmentFlag.AlignHorizontal_Mask
        vertical = flags & Qt.AlignmentFlag.AlignVertical_Mask
        anchor = {
            Qt.AlignmentFlag.AlignRight: "end",
            Qt.AlignmentFlag.AlignLeft: "start",
        }.get(horizontal, "middle")
        vertical_name = {
            Qt.AlignmentFlag.AlignTop: "top",
            Qt.AlignmentFlag.AlignBottom: "bottom",
        }.get(vertical, "middle")
        x = {
            "end": rect.right(),
            "start": rect.left(),
            "middle": rect.center().x(),
        }[anchor]
        y = {
            "top": rect.top(),
            "bottom": rect.bottom(),
            "middle": rect.center().y(),
        }[vertical_name]
        point = mapped(QPointF(x, y))
        svg.text(string, point.x(), point.y(), font, color, anchor, vertical_name)


def _svg_frame(svg, viewbox, rect):
    """The frame of a plot: the edges of a FramedViewBox or the border of a ViewBox."""
    if hasattr(viewbox, "outer_edges"):
        outer = viewbox.outer_edges
        svg.line(
            viewbox.outer_pen if outer["left"] else viewbox.inner_pen,
            rect.topLeft(),
            rect.bottomLeft(),
        )
        svg.line(
            viewbox.outer_pen if outer["top"] else viewbox.inner_pen,
            rect.topLeft(),
            rect.topRight(),
        )
        if outer["right"]:
            svg.line(viewbox.outer_pen, rect.topRight(), rect.bottomRight())
        if outer["bottom"]:
            svg.line(viewbox.outer_pen, rect.bottomLeft(), rect.bottomRight())
    elif viewbox.border is not None:
        stroke = svg.stroke(viewbox.border)
        if stroke:
            svg.add(
                f'<rect x="{rect.x():.2f}" y="{rect.y():.2f}" '
                f'width="{rect.width():.2f}" height="{rect.height():.2f}" {stroke}/>'
            )


def _svg_viewbox(svg, viewbox, to_viewport):
    rect = to_viewport.mapRect(viewbox.sceneBoundingRect())

    # What is outside of the plot is not drawn: the lines are clipped and the markers
    # outside are left out, as not all viewers support clipping
    # (Items that are added with `ignoreBounds` are not in `addedItems`)
    items = sorted(viewbox.childGroup.childItems(), key=lambda item: item.zValue())
    for item in items:
        if not item.isVisible():
            continue
        transform = item.sceneTransform() * to_viewport
        if isinstance(item, pg.PlotDataItem):
            curve = item.curve
            if curve.isVisible():
                _svg_curve(svg, curve, curve.sceneTransform() * to_viewport, rect)
        elif isinstance(item, pg.PlotCurveItem):
            _svg_curve(svg, item, transform, rect)
        elif isinstance(item, pg.ScatterPlotItem):
            _svg_scatter(svg, item, transform, rect)
        elif isinstance(item, pg.InfiniteLine):
            if item.angle in (90, -90):
                x, _ = _map(transform, [item.value()], [0])
                p1, p2 = QPointF(x[0], rect.top()), QPointF(x[0], rect.bottom())
            else:
                _, y = _map(transform, [0], [item.value()])
                p1, p2 = QPointF(rect.left(), y[0]), QPointF(rect.right(), y[0])
            svg.line(item.pen, p1, p2)
        elif isinstance(item, pg.TextItem):
            _svg_text(svg, item, transform)
    _svg_frame(svg, viewbox, rect)


def export_svg(widget, fname):
    """Writes the pyqtgraph plot of the widget to the SVG file.

    The curves, sticks, markers, lines, texts, ticks and frames are written as vector
    graphics. The plots are drawn at their place in the widget."""
    scene = widget.scene()
    to_viewport = widget.viewportTransform()
    size = widget.viewport().size()
    svg = SvgWriter(size.width(), size.height())

    background = widget.backgroundBrush()
    for viewbox in [item for item in scene.items() if isinstance(item, pg.ViewBox)]:
        if viewbox.isVisible():
            _svg_viewbox(svg, viewbox, to_viewport)
    for axis in [item for item in scene.items() if isinstance(item, pg.AxisItem)]:
        if axis.isVisible():
            _svg_axis(svg, axis, to_viewport)

    with open(fname, "w", encoding="utf-8") as file:
        file.write(svg.tostring(background))


def write_figure(widget, fname):
    """Writes the plot to the file: an SVG (vector graphics) or an image like a PNG.

    The format is chosen by the file extension (PNG if there is none). The corners are
    square in both."""
    if fname.lower().endswith(".svg"):
        export_svg(widget, fname)
        return True

    has_extension = bool(os.path.splitext(fname)[1])
    image = widget.grab_clean() if hasattr(widget, "grab_clean") else widget.grab()
    return image.save(fname, None if has_extension else "PNG")


def save_figure_dialog(widget, fname=None):
    """Asks for a file (if no `fname` is given) and saves the plot as a PNG, JPEG or SVG."""
    if fname is None:
        fname = QFileDialog.getSaveFileName(
            widget,
            "Choose file to save the figure to",
            "",
            "PNG Files (*.png);;SVG Files (*.svg);;JPEG Files (*.jpg *.jpeg);;All Files (*)",
        )[0]
    if fname and not write_figure(widget, fname):
        raise OSError(f"The figure could not be saved to '{fname}'.")


def install_save_menu(canvas, reset_view=None):
    """Gives the plot a context menu to save it as a figure.

    If `reset_view` is given, it is called by the entry "Reset View"."""

    def context_menu(event):
        menu = QMenu(canvas)
        reset_action = menu.addAction("Reset View") if reset_view else None
        save_action = menu.addAction("Save Figure")

        action = menu.exec(event.globalPos())
        if action == save_action:
            save_figure_dialog(canvas)
        elif reset_action is not None and action == reset_action:
            reset_view()

    canvas.contextMenuEvent = context_menu


##
## Items
##
class EdgeLabelAxisItem(pg.AxisItem):
    """AxisItem that also draws the labels of ticks at the edges of the plot.

    The labels that are not completely inside of the bounding rectangle are skipped,
    which are the labels at the edges of the range. Therefore it is extended."""

    margin = 60

    def boundingRect(self):
        return super().boundingRect().adjusted(-self.margin, 0, self.margin, 0)


def setup_axis(axis, hide_line=False, outward_ticks=False):
    """Styles an axis with the colors of the theme. Can be called again to update them.

    The line of the axis can be hidden if the frame of the plot is drawn at the same
    place. The ticks can point away from the plot and have all the same opacity."""
    axis.enableAutoSIPrefix(False)
    style = {"tickFont": plot_font()}
    if outward_ticks:
        style.update(tickLength=5, tickAlpha=255)
    axis.setStyle(**style)

    frame = plot_frame_color()
    axis.setPen(pg.mkPen(None) if hide_line else frame)
    axis.setTickPen(frame)
    axis.setTextPen(plot_foreground())


class SpanViewBox(pg.ViewBox):
    """ViewBox that does not zoom or pan by itself, the range is set by the program.

    Dragging with the mouse `button` selects a range, shown as a shaded area. When the
    button is released, the range is passed to `on_span(xmin, xmax)`.

    * `modifier`: keyboard modifier that has to be held to select (by default none)
    * `rectangle`: select a rectangle instead of a range of x, which is passed to
      `on_span(xmin, xmax, ymin, ymax)` if it is at least a few pixels wide and high
    * `limits`: function that returns (xmin, xmax) which the selection is limited to
    * `on_start`: function that is called when the selection starts

    With `native_mouse` the middle button and the wheel still move and zoom the plot as
    usual in pyqtgraph. Dragging the right button or at the axes does nothing.
    """

    def __init__(
        self,
        button=Qt.MouseButton.LeftButton,
        color="#fe6100",
        native_mouse=False,
        modifier=None,
        rectangle=False,
        **kwargs,
    ):
        super().__init__(enableMenu=False, **kwargs)
        self.button = button
        self.native_mouse = native_mouse
        self.modifier = modifier
        self.rectangle = rectangle
        self.on_span = None
        self.on_start = None
        self.limits = None
        self.is_selecting = False
        self.setMouseEnabled(native_mouse, native_mouse)

        brush = pg.mkBrush(to_qcolor(color, 0.25))
        if rectangle:
            self.rubberband = QGraphicsRectItem()
            self.rubberband.setBrush(brush)
            self.rubberband.setPen(pg.mkPen(None))
        else:
            self.rubberband = pg.LinearRegionItem(
                movable=False, brush=brush, pen=pg.mkPen(None)
            )
        self.rubberband.setZValue(200)
        self.rubberband.hide()
        self.addItem(self.rubberband, ignoreBounds=True)

    def mouseDragEvent(self, ev, axis=None):
        if ev.isStart():
            has_modifier = self.modifier is None or bool(ev.modifiers() & self.modifier)
            self.is_selecting = (
                ev.button() == self.button
                and has_modifier
                and self.on_span is not None
                and axis is None
            )

        if not self.is_selecting:
            is_movement = ev.button() == Qt.MouseButton.MiddleButton
            if self.native_mouse and axis is None and is_movement:
                super().mouseDragEvent(ev, axis)
            else:
                ev.ignore()
            return

        ev.accept()
        start, end = self.mapToView(ev.buttonDownPos()), self.mapToView(ev.pos())
        x0, x1 = start.x(), end.x()
        if self.limits is not None:
            xmin, xmax = self.limits()
            x0, x1 = np.clip(x0, xmin, xmax), np.clip(x1, xmin, xmax)
        y0, y1 = start.y(), end.y()

        if ev.isStart() and self.on_start is not None:
            self.on_start()
        if ev.isFinish():
            self.rubberband.hide()
            distance = ev.pos() - ev.buttonDownPos()
            if self.rectangle:
                if abs(distance.x()) >= 6 and abs(distance.y()) >= 6:
                    self.on_span(min(x0, x1), max(x0, x1), min(y0, y1), max(y0, y1))
            else:
                self.on_span(min(x0, x1), max(x0, x1))
        elif self.rectangle:
            self.rubberband.setRect(
                QRectF(
                    QPointF(min(x0, x1), min(y0, y1)), QPointF(max(x0, x1), max(y0, y1))
                )
            )
            self.rubberband.show()
        else:
            self.rubberband.setRegion((min(x0, x1), max(x0, x1)))
            self.rubberband.show()

    def mouseClickEvent(self, ev):
        ev.ignore()

    def wheelEvent(self, ev, axis=None):
        if self.native_mouse:
            super().wheelEvent(ev, axis)
        else:
            ev.ignore()


class FramedViewBox(SpanViewBox):
    """SpanViewBox that draws the frame of a plot that is part of a grid of plots."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.outer_pen = self.inner_pen = pg.mkPen(None)
        self.outer_edges = dict.fromkeys(("left", "top", "right", "bottom"), True)

    def set_frame(self, outer_color=None, inner_color=None, edges=None):
        """Sets the colors of the frame and which of its edges lie at the outside.

        The left and top edge are always drawn, so that the line between two plots is
        only drawn once. The right and bottom edge are only drawn at the outside."""
        if outer_color is not None:
            self.outer_pen = pg.mkPen(outer_color, width=0)
        if inner_color is not None:
            self.inner_pen = pg.mkPen(inner_color, width=0)
        if edges is not None:
            self.outer_edges = edges
        self.update()

    def paint(self, painter, option, widget=None):
        super().paint(painter, option, widget)

        rect = self.rect()
        outer = self.outer_edges
        painter.setPen(self.outer_pen if outer["left"] else self.inner_pen)
        painter.drawLine(rect.topLeft(), rect.bottomLeft())
        painter.setPen(self.outer_pen if outer["top"] else self.inner_pen)
        painter.drawLine(rect.topLeft(), rect.topRight())
        painter.setPen(self.outer_pen)
        if outer["right"]:
            painter.drawLine(rect.topRight(), rect.bottomRight())
        if outer["bottom"]:
            painter.drawLine(rect.bottomLeft(), rect.bottomRight())


class StackedRow:
    """A plot in a column of plots. It remembers the extent of its data."""

    def __init__(self, plot, viewbox):
        self.plot = plot
        self.vb = viewbox
        self.xs, self.ys = [], []
        self.texts = []

    def curve(self, xs, ys, color, width=1, alpha=1.0):
        item = pg.PlotCurveItem(xs, ys, pen=make_pen(color, alpha, width))
        self.plot.addItem(item, ignoreBounds=True)
        self.xs.append(np.asarray(xs))
        self.ys.append(np.asarray(ys))

    def star(self, x, color, y=0):
        item = pg.ScatterPlotItem(
            x=[x], y=[y], size=10, symbol="star", pen=pg.mkPen(None)
        )
        item.setBrush(pg.mkBrush(to_qcolor(color)))
        item.setZValue(100)
        self.plot.addItem(item, ignoreBounds=True)
        self.xs.append(np.array([x]))
        self.ys.append(np.array([y]))

    def vline(self, x, color):
        item = pg.InfiniteLine(x, angle=90, pen=make_pen(color))
        item.setZValue(10)
        self.plot.addItem(item, ignoreBounds=True)
        self.xs.append(np.array([x]))

    def text(self, text, color, font, fraction_x=1, fraction_y=0.95):
        """Text at the top right of the plot, the fractions are relative to its range."""
        item = pg.TextItem(text, color=to_qcolor(color), anchor=(1, 0))
        item.setFont(font)
        item.setZValue(150)
        self.plot.addItem(item, ignoreBounds=True)
        self.texts.append((item, fraction_x, fraction_y))

    def extent(self, values):
        values = np.concatenate(values or [[]])
        return (values.min(), values.max()) if len(values) else None

    def set_ranges(self, xrange=None, yrange=None, ymargin=0.1):
        """Sets the ranges, by default to the extent of the data (y with a margin)."""
        xrange = xrange or self.extent(self.xs)
        yrange = yrange or self.extent(self.ys)
        if xrange is None or yrange is None:
            return

        ymin, ymax = yrange
        if ymin == ymax:
            ymin, ymax = ymin - 1, ymax + 1
        ymin, ymax = ymin - ymargin * (ymax - ymin), ymax + ymargin * (ymax - ymin)
        self.plot.setXRange(*xrange, padding=0)
        self.plot.setYRange(ymin, ymax, padding=0)

        for item, fraction_x, fraction_y in self.texts:
            item.setPos(
                xrange[0] + fraction_x * (xrange[1] - xrange[0]),
                ymin + fraction_y * (ymax - ymin),
            )


def create_stacked_plots(canvas, n_rows, link_x=False, color="#fe6100", tick_density=1):
    """Creates `n_rows` plots in a column without gaps on a GraphicsLayoutWidget.

    The first one is the lowest and the only one with an x axis."""
    canvas.clear()
    frame_color = plot_frame_color()

    rows = []
    for i in range(n_rows):
        viewbox = FramedViewBox(color=color)
        plot = canvas.addPlot(
            row=n_rows - 1 - i,
            col=0,
            viewBox=viewbox,
            axisItems={"bottom": EdgeLabelAxisItem("bottom")},
            enableMenu=False,
        )
        plot.setMenuEnabled(False)
        plot.hideButtons()
        plot.setContentsMargins(0, 0, 0, 0)
        plot.layout.setContentsMargins(0, 0, 0, 0)
        plot.layout.setSpacing(0)
        plot.setMouseEnabled(False, False)
        plot.hideAxis("left")
        if i:
            plot.hideAxis("bottom")
            if link_x:
                plot.setXLink(rows[0].plot)
        else:
            axis = plot.getAxis("bottom")
            setup_axis(axis, hide_line=True, outward_ticks=True)
            axis.setTickDensity(tick_density)

        viewbox.disableAutoRange()
        viewbox.setDefaultPadding(0)
        viewbox.set_frame(
            frame_color,
            plot_inner_frame_color(),
            {"left": True, "top": i == n_rows - 1, "right": True, "bottom": i == 0},
        )
        rows.append(StackedRow(plot, viewbox))
    return rows


def points_in_polygon(xs, ys, polygon_xs, polygon_ys):
    """Boolean mask of the points that are inside of the polygon (even-odd rule)."""
    inside = np.zeros(len(xs), dtype=bool)
    j = len(polygon_xs) - 1
    for i in range(len(polygon_xs)):
        xi, yi, xj, yj = polygon_xs[i], polygon_ys[i], polygon_xs[j], polygon_ys[j]
        crosses = (yi > ys) != (yj > ys)
        with np.errstate(divide="ignore", invalid="ignore"):
            x_crossing = (xj - xi) * (ys - yi) / (yj - yi) + xi
        inside ^= crosses & (xs < x_crossing)
        j = i
    return inside


class LassoViewBox(pg.ViewBox):
    """ViewBox in which the left mouse button zooms to the dragged rectangle.

    * Dragging with the left button zooms to the rectangle. Holding alt moves the plot
      instead (like the middle button), the wheel zooms. Dragging the right button or at
      the axes does nothing.
    * Dragging while holding shift draws a lasso. When the button is released,
      `on_lasso(xs, ys, additive)` is called with the corners of the lasso (`additive`
      is true if control was pressed as well)."""

    def __init__(self, color="#fe6100", **kwargs):
        super().__init__(enableMenu=False, **kwargs)
        self.setMouseMode(pg.ViewBox.RectMode)
        self.on_lasso = None
        self.corners = []
        self.drag_mode = None

        self.lasso_item = pg.PlotCurveItem(
            pen=make_pen(color, style=Qt.PenStyle.DashLine)
        )
        self.lasso_item.setZValue(200)
        self.addItem(self.lasso_item, ignoreBounds=True)

    def mouseDragEvent(self, ev, axis=None):
        # Dragging the right button (also a control-click on a Mac) or dragging at the axes
        # would scale the plot, which is not wanted
        if axis is not None or ev.button() not in (
            Qt.MouseButton.LeftButton,
            Qt.MouseButton.MiddleButton,
        ):
            ev.ignore()
            return

        if ev.isStart():
            modifiers = ev.modifiers()
            is_left = ev.button() == Qt.MouseButton.LeftButton
            if is_left and modifiers & Qt.KeyboardModifier.ShiftModifier:
                self.drag_mode = "lasso" if self.on_lasso is not None else "zoom"
            elif is_left and modifiers & Qt.KeyboardModifier.AltModifier:
                self.drag_mode = "pan"
            else:
                self.drag_mode = "zoom"

        if self.drag_mode == "lasso":
            self.drag_lasso(ev)
        elif self.drag_mode == "pan":
            self.setMouseMode(pg.ViewBox.PanMode)
            try:
                super().mouseDragEvent(ev, axis)
            finally:
                self.setMouseMode(pg.ViewBox.RectMode)
        else:
            if ev.isFinish() and self.is_tiny_drag(ev):
                # A click that moved a little bit must not zoom in to a tiny region
                self.rbScaleBox.hide()
                ev.accept()
                return
            super().mouseDragEvent(ev, axis)

    @staticmethod
    def is_tiny_drag(ev, minimum=6):
        distance = ev.pos() - ev.buttonDownPos()
        return abs(distance.x()) < minimum or abs(distance.y()) < minimum

    def drag_lasso(self, ev):
        ev.accept()
        if ev.isStart():
            self.corners = [self.mapToView(ev.buttonDownPos())]
        self.corners.append(self.mapToView(ev.pos()))
        xs = [corner.x() for corner in self.corners]
        ys = [corner.y() for corner in self.corners]

        if ev.isFinish():
            self.lasso_item.setData([], [])
            if len(xs) >= 3:
                additive = bool(ev.modifiers() & Qt.KeyboardModifier.ControlModifier)
                self.on_lasso(np.array(xs), np.array(ys), additive)
        else:
            self.lasso_item.setData(xs + xs[:1], ys + ys[:1])


class SpectrumItems:
    """The curves and markers of a spectrum (experiment, prediction, assignments).

    The items are reused when the data changes. `apply` has to be called from the GUI thread."""

    def __init__(self, plot, symbol="star", marker_size=10):
        """`symbol` is the name of a symbol of pyqtgraph, for the assignments."""
        self.plot = plot
        self.exp_items, self.cat_items = [], []
        self.lin_item = pg.ScatterPlotItem(
            pxMode=True,
            size=marker_size,
            pen=pg.mkPen(None),
            symbol=symbol,
        )
        self.lin_item.setZValue(100)
        plot.addItem(self.lin_item, ignoreBounds=True)

    def set_curves(self, pool, specs, zvalue=0):
        for i, (xs, ys, connect, color) in enumerate(specs):
            if i == len(pool):
                item = pg.PlotCurveItem(skipFiniteCheck=True)
                item.setZValue(zvalue)
                self.plot.addItem(item, ignoreBounds=True)
                pool.append(item)
            item = pool[i]
            item.setData(xs, ys, connect=connect)
            item.setPen(make_pen(color))
            item.show()
        for item in pool[len(specs) :]:
            item.hide()

    def set_markers(self, xs, ys, colors):
        brushes = {color: pg.mkBrush(to_qcolor(color)) for color in set(colors)}
        self.lin_item.setData(
            x=xs, y=ys, brush=[brushes[color] for color in colors], pen=pg.mkPen(None)
        )

    def apply(self, spec):
        self.set_curves(self.exp_items, spec["exp"])
        self.set_curves(self.cat_items, spec["cat"], zvalue=1)
        self.set_markers(*spec["lin"])


class ScatterPlotWidget(RoundedPlotWidget):
    """Scatter plot with a tooltip for the hovered points that can be zoomed and moved.

    * `hovered(indices)` is emitted when the mouse moves over the plot, `indices` are the
      indices of the points under the cursor (empty if there are none)
    * `clicked(indices, button)` is emitted when the plot is clicked
    * `selection_changed(indices)` is emitted when the selection changed. Points are
      selected by drawing a lasso with shift pressed (additionally with control to add
      to the selection). A click on an empty spot clears the selection.

    Dragging a rectangle zooms to it (alt moves the plot instead, see `LassoViewBox`),
    the context menu resets the view to show all points.

    The tooltip is shown with `show_tooltip`. `set_points` has to be called from the
    GUI thread. The context menu lists the actions of `context_actions(hovered, selected)`,
    which returns a list of (text, function), and the saving of the figure."""

    hovered = pyqtSignal(object)
    clicked = pyqtSignal(object, object)
    selection_changed = pyqtSignal(object)

    def __init__(self, marker_size=6):
        self.vb = LassoViewBox()
        super().__init__(viewBox=self.vb, background=plot_background())
        self.marker_size = marker_size
        self.ymargin = 0.1
        self.xs = self.ys = np.array([])
        self.selected = np.zeros(0, dtype=bool)
        self.context_actions = None
        self.vb.on_lasso = self.select_in_lasso

        self.plot = self.getPlotItem()
        self.plot.setMenuEnabled(False)
        self.plot.hideButtons()
        for name in ("bottom", "left"):
            setup_axis(self.plot.getAxis(name))

        self.scatter = pg.ScatterPlotItem(
            size=marker_size, symbol="o", pen=pg.mkPen(None)
        )
        self.plot.addItem(self.scatter)

        self.selection_item = pg.ScatterPlotItem(
            size=marker_size + 6,
            symbol="o",
            pen=pg.mkPen("#ff0000", width=2),
            brush=pg.mkBrush(None),
        )
        self.selection_item.setZValue(150)
        self.plot.addItem(self.selection_item, ignoreBounds=True)

        self.tooltip = pg.TextItem(
            color="#000000",
            anchor=(0.5, 1.0),
            fill=pg.mkBrush("#ffffff"),
            border=pg.mkPen("#000000"),
        )
        self.tooltip.setFont(plot_font())
        self.tooltip.setZValue(200)
        self.tooltip.hide()
        self.plot.addItem(self.tooltip, ignoreBounds=True)

        self.scene().sigMouseMoved.connect(self.on_mouse_moved)
        self.scene().sigMouseClicked.connect(self.on_mouse_clicked)
        theme.changed.connect(self.apply_theme)
        self.apply_theme()

    def apply_theme(self):
        self.setBackground(plot_background())
        self.vb.setBorder(pg.mkPen(plot_frame_color()))
        for name in ("bottom", "left"):
            setup_axis(self.plot.getAxis(name))

    def set_points(self, xs, ys, colors, autoscale=False, ymargin=0.1):
        self.xs, self.ys = np.asarray(xs, dtype=float), np.asarray(ys, dtype=float)
        brushes = {color: pg.mkBrush(to_qcolor(color)) for color in set(colors)}
        self.scatter.setData(
            x=self.xs,
            y=self.ys,
            brush=[brushes[color] for color in colors],
            pen=pg.mkPen(None),
        )
        self.tooltip.hide()
        self.set_selection(np.zeros(len(self.xs), dtype=bool))

        self.ymargin = ymargin
        if autoscale:
            self.reset_view()

    def reset_view(self):
        """Shows all points."""
        if not len(self.xs):
            return

        xmin, xmax = self.xs.min(), self.xs.max()
        ymin, ymax = self.ys.min(), self.ys.max()
        if xmin == xmax:
            xmin, xmax = xmin - 1, xmax + 1
        if ymin == ymax:
            ymin, ymax = ymin - 1, ymax + 1
        margin = self.ymargin
        self.plot.setXRange(
            xmin - 0.02 * (xmax - xmin), xmax + 0.02 * (xmax - xmin), padding=0
        )
        self.plot.setYRange(
            ymin - margin * (ymax - ymin), ymax + margin * (ymax - ymin), padding=0
        )

    # Selection
    def set_selection(self, selected):
        self.selected = selected
        self.selection_item.setData(x=self.xs[selected], y=self.ys[selected])
        self.selection_changed.emit(self.selected_indices())

    def selected_indices(self):
        return np.flatnonzero(self.selected)

    def clear_selection(self):
        if self.selected.any():
            self.set_selection(np.zeros(len(self.xs), dtype=bool))

    def select_in_lasso(self, lasso_xs, lasso_ys, additive):
        inside = points_in_polygon(self.xs, self.ys, lasso_xs, lasso_ys)
        self.set_selection(self.selected | inside if additive else inside)

    # Tooltip and mouse
    def show_tooltip(self, x, y, text):
        self.tooltip.setText(text)
        self.tooltip.setPos(x, y)
        self.tooltip.show()

    def hide_tooltip(self):
        self.tooltip.hide()

    def points_at(self, scene_pos):
        """Indices of the points under the position in the scene."""
        if not self.vb.sceneBoundingRect().contains(scene_pos) or not len(self.xs):
            return np.array([], dtype=int)

        point = self.vb.mapSceneToView(scene_pos)
        pixel_width, pixel_height = self.vb.viewPixelSize()
        radius = self.marker_size / 2 + 2
        near = (np.abs(self.xs - point.x()) < radius * pixel_width) & (
            np.abs(self.ys - point.y()) < radius * pixel_height
        )
        return np.flatnonzero(near)

    def on_mouse_moved(self, scene_pos):
        indices = self.points_at(scene_pos)
        if not len(indices):
            self.hide_tooltip()
        self.hovered.emit(indices)

    def on_mouse_clicked(self, event):
        indices = self.points_at(event.scenePos())
        if event.button() == Qt.MouseButton.LeftButton and not len(indices):
            self.clear_selection()
        self.clicked.emit(indices, event.button())

    def contextMenuEvent(self, event):
        hovered = self.points_at(self.mapToScene(event.pos()))
        actions = []
        if self.context_actions is not None:
            actions = self.context_actions(hovered, self.selected_indices())

        menu = QMenu(self)
        functions = {menu.addAction(text): function for text, function in actions}
        if actions:
            menu.addSeparator()
        reset_action = menu.addAction("Reset View")
        save_action = menu.addAction("Save Figure")

        action = menu.exec(event.globalPos())
        if action == reset_action:
            self.reset_view()
        elif action == save_action:
            save_figure_dialog(self)
        elif action in functions:
            functions[action]()
