import cairo 
import re
import numpy as np


class drawer():
    # The class utilizies drawing using cairo library. 
    # Here, different graphic shapes can be defined within separate methods.

    def __init__(self, out_path = None, width = 1980, height = 1224, bg_color = (1,1,1)):
        
        # input handling
        
        if (out_path == None):
            raise NameError("There must be out_path specified!")
        else:
            self.out_path = out_path
        
        # creating cairo surface for drawing
        
        path_feature = self.path_features(out_path)
        
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
        self.ctx.rectangle(0,0,width,height)
        self.ctx.fill()


    def path_features(self,filePath):
        # this function fetches features from provided path
        matching_pattern = "((?:(?:.*\\/)*(.*)\\/))*(?:(.*)\\.(.*))*"
        feature_match = re.match(matching_pattern, filePath)
        out_dict = {"parent_dir": feature_match.group(2),
                    "file_name": feature_match.group(3),
                    "extension": feature_match.group(4),
                    "path": feature_match.group(1)}
        return out_dict


    def draw_text(self, text = None, x = None, y = None, font_size = 28, rotate = 0, color_rgb = (0,0,0), h_align = "center", v_align = "center"):
        # This function draws text according to the specifications to the cairo context
        
        if (text == None):
            return
        
        self.ctx.set_source_rgb(color_rgb[0],color_rgb[1],color_rgb[2])
        self.ctx.set_font_size(font_size)
        
        (_, _, width, height, _, _) = self.ctx.text_extents(text) 

        if h_align == "center":
            x = x - width/2 
        elif h_align == "right":
            x = x - width
        elif h_align == "left":
            x = x

        if v_align == "center":
            y = y + height/2
        elif v_align == "bottom":
            y = y + height
        elif v_align == "top":
            y = y

        self.ctx.move_to(x,y)
        self.ctx.save()
        self.ctx.rotate(rotate * np.pi/180)
        self.ctx.show_text(text)
        self.ctx.restore()
        #self.ctx.rotate(-rotate * np.pi/180)    

    def draw_line(self, x0 = None, y0 = None, width = None, end_shape = ("both","arrow"), color = (0,0,0), line_width = 4, rotate = 0, alpha = 1):
        # Draw arrow line
        
        ctx = self.ctx
        ctx.set_source_rgba(color[0],color[1],color[2], alpha)
        ctx.set_line_width(line_width)
        ctx.set_line_cap(cairo.LINE_CAP_ROUND)
        
        ctx.save()
        ctx.translate(x0, y0)
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
            if (end_shape[0] in ("both","left")):
                ctx.move_to(rel_x + arrow_x, rel_y - arrow_y)
                ctx.line_to(rel_x, rel_y)
                ctx.line_to(rel_x + arrow_x, rel_y + arrow_y)
            
            if (end_shape[0] in ("both","right")):
                ctx.move_to(x - arrow_x, rel_y - arrow_y)
                ctx.line_to(x, rel_y)
                ctx.line_to(x - arrow_x, rel_y + arrow_y)
        
        elif (end_shape[1] == "line"):
            line_length = 8                
            
            if (end_shape[0] in ("both","left")):
                ctx.move_to(rel_x, rel_y - line_length/2)
                ctx.line_to(rel_x, rel_y + line_length/2)
            
            if (end_shape[0] in ("both","right")):
                ctx.move_to(x, rel_y - line_length/2)
                ctx.line_to(x, rel_y + line_length/2)
        
        ctx.stroke()
        ctx.restore()
        
    def draw_rectangle(self, x = None, y = None, width = None, height = None, round_aspect = None, line_width = 4, color_line = (0,0,0), color_fill = (0,0,0), rotate = 0, alpha = 1):
        ctx = self.ctx
        
        ctx.set_line_width(line_width)
        ctx.save()
        ctx.translate(x, y)
        
        if (rotate != 0):
            ctx.rotate(rotate*np.pi/180)
        
        # Fix 0 division 
        if round_aspect == 0:
            round_aspect = None
        
        if round_aspect != None:
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
            ctx.rectangle(0,0,width,height)
            
        ctx.set_source_rgba(color_fill[0],color_fill[1],color_fill[2],alpha)
        ctx.fill_preserve()    
        ctx.set_source_rgba(color_line[0],color_line[1],color_line[2],alpha)
        ctx.stroke()
        ctx.restore()
        
    def return_differences(self, vector):
        ## This function gives indexes of non-zero differences
        ## input: ndarray num vector
        ## output: nddarray of position of differences
        differences = np.argwhere((np.diff(vector)!=0).astype(int)).flatten()    
        differences = np.append(0,differences)
        differences = np.append(differences,len(vector))
        
        return differences
        
    def draw_signal(self, signal = None, x = None, y = None, width = None, height = None, y_max = None, mode = "normal", color_fill = (0.3,0.3,0.3), color_line = (0,0,0), line_width = 4, rotate = 0):
        ctx = self.ctx
        ctx.set_line_cap(cairo.LINE_CAP_ROUND)
        ctx.set_line_width(line_width)
        
        ctx.save()
        ctx.translate(x, y + height/2)
        ctx.rotate(rotate * np.pi/180)
        
        sig_length = len(signal)
        
        steps = np.linspace(len(signal)/width, width, len(signal))
        
        rel_x = 0
        rel_y = 0
        
        ctx.new_path()
        ctx.move_to(rel_x, rel_y) 
        
        if y_max == None:
            y_max = max(signal)
    
        calculate_pos_y = lambda x: x / y_max * height
        val_y = calculate_pos_y(signal[0])
        
        ctx.line_to(rel_x, rel_y - val_y)
        
        if (mode == "normal"):
            for index in range(1,len(signal)):
                val_y = calculate_pos_y(signal[index])
                ctx.line_to(rel_x + steps[index],rel_y - val_y)
        
        elif (mode == "segment"):
            vector = self.return_differences(signal)
            
            for index_diff, index_sig in enumerate(vector[:-1]):
                
                ctx.line_to(rel_x + steps[index_sig], rel_y - val_y)
                val_y = calculate_pos_y(signal[index_sig + 1])
                ctx.line_to(rel_x + steps[index_sig], rel_y - val_y)
            
            ctx.line_to(rel_x + width, rel_y - val_y)
        
        ctx.line_to(rel_x + width, rel_y)
        
        ctx.set_line_width(1)
        ctx.close_path()
        ctx.set_source_rgba(color_fill[0],color_fill[1],color_fill[2],1)
        ctx.set_source_rgba(0.5,0,1,1)
        ctx.fill_preserve()
        ctx.set_source_rgba(color_line[0],color_line[1],color_line[2],1)
        ctx.stroke()
        ctx.restore()
        
            
    def finish(self):
        ext = self.path_features(self.out_path)["extension"]
        surface = self.surface
        ctx = self.ctx

        surface.show_page()
        
        # Handling for ImageSurface is needed
        if (ext == "png"):
            surface.write_to_png(self.out_path)

        surface.finish()
        
        del ctx