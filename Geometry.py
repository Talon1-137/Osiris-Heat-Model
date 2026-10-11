# The DXF controls the hot-wall radius. Channel dimensions below are in mm.
# Flow areas, hydraulic diameter, and heated area are returned in SI units.

import math


def channel_geometry(x, contour, inputs):
    throat = contour.throat()
    radius = contour.radius_at(x)

    if x <= throat["x"]:
        end_radius = contour.radius_at(contour.x_min)
        end_width = inputs["chamber_channel_width"]
    else:
        end_radius = contour.radius_at(contour.x_max)
        end_width = inputs["exit_channel_width"]

    fraction = (radius - throat["radius"]) / (end_radius - throat["radius"])
    width = (
        inputs["throat_channel_width"]
        + fraction * (end_width - inputs["throat_channel_width"])
    ) / 1000
    height = inputs["channel_height"] / 1000
    hot_perimeter = 2 * math.pi * radius / 1000
    channel_pitch = (
        2 * math.pi * (radius + inputs["wall_thickness"])
        / 1000 / inputs["number_of_channels"]
    )
    if width >= channel_pitch:
        raise ValueError(f"Cooling channels overlap at x={x:.3f} mm")

    return {
        "width": width,
        "height": height,
        "rib_width": channel_pitch - width,
        "hot_wall_pitch": hot_perimeter / inputs["number_of_channels"],
        "flow_area": width * height,
        "hydraulic_diameter": 4 * width * height / (2 * (width + height)),
    }


def wall_segment(x_start, x_end, contour):
    radius_start = contour.radius_at(x_start) / 1000
    radius_end = contour.radius_at(x_end) / 1000
    length = math.hypot(
        (x_end - x_start) / 1000, radius_end - radius_start
    )
    # Lateral area of a conical frustum, including the contour slope.
    return math.pi * (radius_start + radius_end) * length


def axial_stations(contour, count):
    if not isinstance(count, int) or isinstance(count, bool) or count < 1:
        raise ValueError("station_count must be a positive integer")

    # Include the throat as a cell boundary so no segment crosses the Mach branch.
    edges = [
        contour.x_min + i * (contour.x_max - contour.x_min) / count
        for i in range(count + 1)
    ]
    edges[0], edges[-1] = contour.x_min, contour.x_max
    throat_x = contour.throat()["x"]
    if not any(math.isclose(x, throat_x, abs_tol=1e-9) for x in edges):
        edges.append(throat_x)
    edges.sort()

    return [
        {
            "x": (start + end) / 2,
            "heated_area": wall_segment(start, end, contour),
        }
        for start, end in zip(edges, edges[1:])
    ]
