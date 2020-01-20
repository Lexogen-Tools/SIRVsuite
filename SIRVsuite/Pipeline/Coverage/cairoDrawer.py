import cairo
import numpy as np
import os
import logging
from ..helper import *

log = logging.getLogger(__name__.split(".")[-1])

class CairoDrawer():
    """
    The class utilizies drawing using cairo library.
    Here, different graphic shapes can be defined within separate methods.
    """

    def __init__(self, out_path=None, width=2000, height=1324, bg_color=(1, 1, 1)):
        """
        During the initilization, define surface format, dimensions and background color (manifested by a rectangle of the surface size)
        """
        # input handling

        if (out_path is None):
            raise NameError("There must be out_path specified!")
        else:
            self.out_path = out_path

        path = path_features(out_path)["path"]
        if not os.path.exists(path):
            os.makedirs(path)

        # creating cairo surface for drawing

        path_feature = path_features(out_path)

        if path_feature["extension"] == "pdf":
            self.surface = cairo.PDFSurface(out_path, width, height)
        elif path_feature["extension"] == "svg":
            self.surface = cairo.SVGSurface(out_path, width, height)
        elif path_feature["extension"] == "ps":
            self.surface = cairo.PSSurface(out_path, width, height)

        # different syntax for ImageSurface, out_path provided before finish
        elif path_feature["extension"] == "png":
            self.surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, width, height)

        else:
            raise NameError("Unsupported graphics format.. Try to use .svg, .pdf, .png or .ps extension.")

        # creating background rectangle to avoid transparent background in the final result

        self.ctx = cairo.Context(self.surface)
        self.ctx.set_source_rgb(bg_color[0], bg_color[1], bg_color[2])
        self.ctx.rectangle(0, 0, width, height)
        self.ctx.fill()

    def get_text_size(self, text=None, fontsize=None):
        """
        A helper method to return text dimensions
        """
        self.ctx.save()
        self.ctx.set_font_size(fontsize)
        obj = self.ctx.text_extents(text)
        self.ctx.restore()
        return obj.width, obj.height

    def __merge_close_coordinates__(self, thre=50, starts=None, ends=None):
        """
        A helper method to prevent overlap of mark triangle symbols for exceeding coverage in coverage plot
        """
        if (len(starts) == 0 or len(ends) == 0):
            return starts, ends

        idx = 0

        starts_merged = []
        ends_merged = []

        while idx < len(ends):
            starts_merged.append(starts[idx])
            pos_thre = starts[idx] + thre
            end_idx_lower = np.argwhere(ends <= pos_thre).flatten()
            start_idx_lower = np.argwhere(starts <= pos_thre).flatten()

            if len(end_idx_lower) == 0 and len(start_idx_lower) == 0:
                end_idx = idx
            else:
                end_idx = max(np.append(end_idx_lower, start_idx_lower))

            ends_merged.append(ends[end_idx])

            idx = end_idx + 1

        return starts_merged, ends_merged

    def draw_text(self, text=None, x=None, y=None, font_size=28, rotate=0, color_rgb=(0, 0, 0), alpha=1, h_align="center", v_align="center", rotation_point=""):
        """
        A method to draw text onto the surface
        """

        self.ctx.save()

        if (text is None):
            return

        self.ctx.set_source_rgba(color_rgb[0], color_rgb[1], color_rgb[2], alpha)
        self.ctx.set_font_size(font_size)

        txt_obj = self.ctx.text_extents(text)

        if h_align == "center":
            dx = -txt_obj.width/2
        elif h_align == "right":
            dx = -txt_obj.width
        elif h_align == "left":
            dx = 0

        if v_align == "center":
            dy = txt_obj.height/2 - (txt_obj.height + txt_obj.y_bearing)
        elif v_align == "bottom":
            dy = txt_obj.height/2 + (txt_obj.height + txt_obj.y_bearing)
        elif v_align == "top":
            dy = 0

        self.ctx.translate(x, y)

        if rotate not in [0, 180, 360]:
            x_shift = -dy
        else:
            x_shift = 0

        if rotation_point == "center":
            self.ctx.rotate(rotate * np.pi/180)
            self.ctx.translate(dx, dy)
        else:
            self.ctx.translate(dx+x_shift, dy)
            self.ctx.rotate(rotate * np.pi/180)

        self.ctx.move_to(0, 0)

        self.ctx.show_text(text)
        self.ctx.restore()

    def draw_line(self, x=None, y=None, width=None, end_shape=("", ""), color=(0, 0, 0), line_width=4, rotate=0, alpha=1):
        """
        A method to draw line onto the surface
        """
        # Draw arrow line

        ctx = self.ctx
        ctx.set_source_rgba(color[0], color[1], color[2], alpha)
        ctx.set_line_width(line_width)
        ctx.set_line_cap(cairo.LINE_CAP_ROUND)

        ctx.save()
        ctx.translate(x, y)
        ctx.rotate(rotate * np.pi/180.0)

        rel_x = 0
        rel_y = 0
        x = rel_x + width

        ctx.move_to(rel_x, rel_y)
        ctx.line_to(x, rel_y)

        if (end_shape[1] == "arrow"):
            # predefined arrow size
            arrow_y = 6
            arrow_x = 8

            # draw arrow
            if (end_shape[0] in ("both", "left")):
                ctx.move_to(rel_x + arrow_x, rel_y - arrow_y)
                ctx.line_to(rel_x, rel_y)
                ctx.line_to(rel_x + arrow_x, rel_y + arrow_y)

            if (end_shape[0] in ("both", "right")):
                ctx.move_to(x - arrow_x, rel_y - arrow_y)
                ctx.line_to(x, rel_y)
                ctx.line_to(x - arrow_x, rel_y + arrow_y)

        elif (end_shape[1] == "line"):
            line_length = 8

            if (end_shape[0] in ("both", "left")):
                ctx.move_to(rel_x, rel_y - line_length/2)
                ctx.line_to(rel_x, rel_y + line_length/2)

            if (end_shape[0] in ("both", "right")):
                ctx.move_to(x, rel_y - line_length/2)
                ctx.line_to(x, rel_y + line_length/2)

        ctx.stroke()
        ctx.restore()

    def draw_rectangle(self, x=None, y=None, width=None, height=None, round_aspect=None, line_width=0, color_line=(0, 0, 0), color_fill=(0, 0, 0), rotate=0, alpha=1):
        """
        A method to draw rectangle onto the surface
        """
        ctx = self.ctx

        ctx.set_line_width(line_width)
        ctx.save()
        ctx.translate(x, y)

        if (rotate != 0):
            ctx.rotate(rotate*np.pi/180)

        # Fix 0 division
        if round_aspect == 0:
            round_aspect = None

        if round_aspect is not None:
            corner_radius = height / 10.0
            radius = corner_radius / float(round_aspect)
            degrees = np.pi / 180.0
            relative_x = 0
            relative_y = 0

            ctx.new_sub_path()
            ctx.arc(relative_x + width - radius, relative_y + radius, radius, -90 * degrees, 0 * degrees)
            ctx.arc(relative_x + width - radius, relative_y + height - radius, radius, 0 * degrees, 90 * degrees)
            ctx.arc(relative_x + radius, relative_y + height - radius, radius, 90 * degrees, 180 * degrees)
            ctx.arc(relative_x + radius, relative_y + radius, radius, 180 * degrees, 270 * degrees)
            ctx.close_path()
        else:
            ctx.rectangle(0, 0, width, height)

        ctx.set_source_rgba(color_fill[0], color_fill[1], color_fill[2], alpha)
        ctx.fill_preserve()
        ctx.set_source_rgba(color_line[0], color_line[1], color_line[2], 1)
        ctx.stroke()
        ctx.restore()

    def draw_triangle(self, x=None, y=None, width=None, height=None, rotate=0, color_fill=None, color_line=(0, 0, 0), alpha=1, line_width=4):
        """
        A method to draw a triangle onto the surface
        """
        ctx = self.ctx

        ctx.save()
        ctx.set_line_width(line_width)
        ctx.set_line_cap(cairo.LINE_CAP_ROUND)
        ctx.translate(x, y)
        ctx.rotate(rotate * np.pi/180)

        rel_x = 0
        rel_y = 0

        ctx.move_to(rel_x - width/2, rel_y + height/2)
        ctx.line_to(rel_x + width/2, rel_y + height/2)
        ctx.line_to(rel_x, rel_y - height/2)
        ctx.close_path()

        if (color_fill is not None):
            ctx.set_source_rgba(color_fill[0], color_fill[1], color_fill[2], alpha)
            ctx.fill_preserve()

        ctx.set_source_rgba(color_line[0], color_line[1], color_line[2], 1)
        ctx.stroke()
        ctx.restore()

    def return_differences(self, vector, boundary_fill=False, mode='all'):
        """
        This method gives indexes of non-zero, positive or negative differences from inserted list or numpy array
        """

        if (not boundary_fill):
            vector = np.append(0, vector)
            vector = np.append(vector, 0)

        if (mode == 'all'):
            differences = np.argwhere((np.diff(vector) != 0).astype(int)).flatten()
        elif (mode == 'positive'):
            differences = np.argwhere((np.diff(vector) > 0).astype(int)).flatten()
        elif (mode == 'negative'):
            differences = np.argwhere((np.diff(vector) < 0).astype(int)).flatten()

        if (boundary_fill):
            differences = np.append(0, differences)
            differences = np.append(differences, len(vector))

        return differences

    def draw_signal(self, signal=None, x=None, y=None, width=None, height=None, y_max=None, mode="normal", color_fill=(0.3, 0.3, 0.3), color_line=(0, 0, 0), line_width=4, rotate=0, alpha_fill=1, alpha_line=1, upside_down=False, scale_factor=1):
        """
        A method to draw 1D signal(numpy array, list) within a rectangular area specified by width, height, x and y position parameters. 
        Drawing method can be chosen from normal and segment mode.
        
        Segment mode:
            Searches for invariant parts of a signal, any change of difference is displayed on the same x position (creates star-like plot). This was designed for theoretical coverage plot.
        Normal mode:    
            Plots consecutively all values from signal nevertheless the differences. 

        Parameters:
            y_max - defines maximal value, which is plotted at maximal height of the rectangular plot.
            scale_factor - scales the values which exceeds the y_max value
        """
        ctx = self.ctx

        ctx.save()
        ctx.set_line_cap(cairo.LINE_CAP_ROUND)
        ctx.set_line_width(line_width)
        ctx.translate(x, y)
        ctx.rotate(rotate * np.pi/180)

        if (upside_down):
            ctx.scale(1, -1)

        signal = np.array(signal)
        steps = np.linspace(0, width, len(signal))
        steps = np.append(steps, steps[-1])

        rel_x = 0
        rel_y = 0

        ctx.new_path()
        can_close_path = True

        if y_max == 0:
            y_max = 1

        if y_max is not None:
            exceeding_part = np.array(signal > y_max).astype(int)

            exceed_starts = self.return_differences(exceeding_part, mode='positive')
            exceed_ends = self.return_differences(exceeding_part, mode='negative')

            starts_merged, ends_merged = self.__merge_close_coordinates__(starts=exceed_starts, ends=exceed_ends, thre=100)

            max_exceed = []
            for i in range(len(starts_merged)):
                max_exceed.append(max(signal[starts_merged[i]:ends_merged[i]])*scale_factor)

            arrow_width = 15
            arrow_height = 10
            gap_arrow_y = arrow_height
            
            signal[signal >= y_max] = y_max

        if y_max is not None:
            def calculate_pos_y(x): return x / y_max * height
        else:
            def calculate_pos_y(x): return x / max(signal) * height

        val_y = calculate_pos_y(signal[0])

        if (mode == "normal"):

            ctx.move_to(rel_x, rel_y)
            ctx.line_to(rel_x, rel_y - val_y)

            for index in range(1, len(signal)):
                val_y = calculate_pos_y(signal[index])
                ctx.line_to(rel_x + steps[index], rel_y - val_y)

            ctx.line_to(rel_x + width, rel_y)

        elif (mode == "segment"):
            vector = self.return_differences(signal, boundary_fill=False)

            if len(vector) != 0:

                ctx.move_to(rel_x + steps[vector[0]], rel_y)
                ctx.line_to(rel_x + steps[vector[0]], rel_y - val_y)

                for index_diff, index_sig in enumerate(vector[0:-1]):

                    ctx.line_to(rel_x + steps[index_sig], rel_y - val_y)
                    val_y = calculate_pos_y(signal[index_sig])
                    ctx.line_to(rel_x + steps[index_sig], rel_y - val_y)

                ctx.line_to(rel_x + steps[vector[-1]], rel_y - val_y)
                ctx.line_to(rel_x + steps[vector[-1]], rel_y)
            else:
                can_close_path = False

        if (can_close_path):
            ctx.close_path()
            ctx.set_source_rgba(color_fill[0], color_fill[1], color_fill[2], alpha_fill)
            ctx.fill_preserve()

        ctx.set_source_rgba(color_line[0], color_line[1], color_line[2], alpha_line)
        ctx.stroke()
        ctx.restore()

        if (y_max is not None):

            if (not upside_down):
                triangle_y = y - height - gap_arrow_y
                line_y = y - height
                rotate = 0
                factor = -1
            else:
                triangle_y = y + height + gap_arrow_y
                line_y = y + height
                factor = 1
                rotate = 180

            for i in range(len(exceed_starts)):
                self.draw_line(x=x + steps[exceed_starts[i]], y=line_y, width=(steps[exceed_ends[i]-1]-steps[exceed_starts[i]]), color=(1, 0, 0), line_width=2)

            for i in range(len(starts_merged)):
                self.draw_triangle(x=x + steps[starts_merged[i]] + (steps[ends_merged[i]-1]-steps[starts_merged[i]])/2,
                                   y=triangle_y,
                                   width=arrow_width,
                                   height=10,
                                   color_line=(0.5, 0, 0),
                                   line_width=2,
                                   alpha=1,
                                   rotate=rotate)

                self.draw_text(
                    text=str(int(max_exceed[i])), x=x + steps[starts_merged[i]] + (steps[ends_merged[i]]-steps[starts_merged[i]])/2,
                    y=triangle_y + factor*arrow_width, v_align="center",
                    h_align="center", color_rgb=(0.5, 0, 0), font_size=12
                    )

    def draw_table(self, table_dict=None, x=None, y=None, width=None, height=None):
        """
        A method to draw table of attribute-value pairs separated by colon in the middle from a dictionary input.
        """

        num_lines = len(table_dict)
        atributes = table_dict.keys()
        values = table_dict.values()

        row_height = height / num_lines
        font_size = height / num_lines * .7
        text_y = y + row_height/2

        text_width, text_height = self.get_text_size(text=list(atributes)[0], fontsize=font_size)

        for i in range(num_lines):
            self.draw_text(text=list(atributes)[i]+":", x=x+width/2, y=text_y, font_size=font_size, h_align="right", v_align="center")
            self.draw_text(text=list(values)[i], x=x+width/2 + 20, y=text_y, font_size=font_size, h_align="left", v_align="center")
            text_y += row_height

    def finish(self):
        """
        A method to render the surface and save to the specified output.
        """
        ext = path_features(self.out_path)["extension"]
        surface = self.surface
        ctx = self.ctx

        surface.show_page()

        # Handling for ImageSurface is needed
        if (ext == "png"):
            surface.write_to_png(self.out_path)

        surface.finish()

        del ctx
