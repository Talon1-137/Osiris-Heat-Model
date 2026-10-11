# Reads the RPA DXF profile and calculates axisymmetric flow area along x.
# The DXF coordinates and the x argument are in mm; returned areas are in m^2.

import math
from pathlib import Path

import ezdxf


class Contour:
    def __init__(self, segments):
        self.segments = segments
        if not segments:
            raise ValueError("The DXF did not contain an upper chamber/nozzle profile")

        self.x_min = min(min(a[0], b[0]) for a, b in segments)
        self.x_max = max(max(a[0], b[0]) for a, b in segments)

    def radius_at(self, x):
        if x < self.x_min or x > self.x_max:
            raise ValueError(f"x={x} mm is outside the contour ({self.x_min} to {self.x_max} mm)")

        radii = []
        for (x0, y0), (x1, y1) in self.segments:
            low_x = min(x0, x1)
            high_x = max(x0, x1)
            if x < low_x or x > high_x:
                continue

            if math.isclose(x0, x1, abs_tol=1e-12):
                if math.isclose(x, x0, abs_tol=1e-9):
                    radii.extend((y0, y1))
                continue

            fraction = (x - x0) / (x1 - x0)
            radius = y0 + fraction * (y1 - y0)
            if radius > 0:
                radii.append(radius)

        if not radii:
            raise ValueError(f"No contour profile found at x={x} mm")

        # Use the outermost positive-y wall in case adjoining DXF entities overlap.
        return max(radii)

    def area_at(self, x):
        radius_m = self.radius_at(x) / 1000
        return math.pi * radius_m**2

    def throat(self):
        points = [point for segment in self.segments for point in segment]
        x, radius = min(points, key=lambda point: point[1])
        return {
            "x": x,
            "radius": radius,
            "area": math.pi * (radius / 1000)**2,
        }


def load_contour(path=None):
    if path is None:
        path = Path(__file__).with_name("contour.dxf")
    doc = ezdxf.readfile(path)
    if doc.header.get("$INSUNITS") != 4:
        raise ValueError("The RPA contour DXF must use millimeter units ($INSUNITS = 4)")

    segments = []
    for entity in doc.modelspace():
        kind = entity.dxftype()
        points = []

        if kind == "LINE":
            start = entity.dxf.start
            end = entity.dxf.end
            points = [(start.x, start.y), (end.x, end.y)]
        elif kind == "ARC":
            center = entity.dxf.center
            radius = entity.dxf.radius
            start_angle = entity.dxf.start_angle
            end_angle = entity.dxf.end_angle
            while end_angle < start_angle:
                end_angle += 360
            segment_count = max(2, math.ceil((end_angle - start_angle) / 2))
            for index in range(segment_count + 1):
                angle = math.radians(
                    start_angle + (end_angle - start_angle) * index / segment_count
                )
                points.append((
                    center.x + radius * math.cos(angle),
                    center.y + radius * math.sin(angle),
                ))
        elif kind == "SPLINE":
            points = [(point.x, point.y) for point in entity.flattening(distance=0.01)]

        # The RPA file contains the upper and lower mirrored walls. Keep only
        # positive-y geometry and ignore the centerline and dimensions.
        for start, end in zip(points, points[1:]):
            if start[1] > 0 and end[1] > 0:
                segments.append((start, end))

    return Contour(segments)
