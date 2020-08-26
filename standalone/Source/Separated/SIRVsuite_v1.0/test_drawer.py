from Pipeline.Coverage.CairoDrawer import *

def get_continous_coverage_ends(starts, ends):
    #transcripts = coord_dict.keys()
    #exon_starts = np.array([coord_dict[i]["exon_starts"] for i in transcripts]).flatten()
    #exon_ends = np.array([coord_dict[i]["exon_ends"] for i in transcripts]).flatten()
    
    sorted_idx = np.argsort(starts)

    contig_start = [starts[sorted_idx[0]]]
    contig_end = []

    for i in range(0,len(sorted_idx)-1):

        prev_idx = ends[sorted_idx[i]]
        next_idx = starts[sorted_idx[i+1]]

        if (next_idx > prev_idx):
            contig_end.append(prev_idx)
            contig_start.append(next_idx)

        
    contig_end.append(ends[sorted_idx[-1]])


    return contig_start, contig_end

signal_segment = [5,5,5,5,10,10,10,10,10,10,10,10,0,0,0,0,0,0,2,2,2,2,4,4,6,6,0,0,0,0,0,4,4,4,8,8,8,8,15,15,15,15,15,15,15,15,15]
signal_normal = [45,55,53,10,0,0,10,36,40,100,50,20,0,0,0,0,0,0,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,14,13,12,11,10,9,8,7,6,5,4,59,59,59,59]
t_name = ["SIRV1","SIRV2"]
t_name = t_name[0]
transcript1 = np.array([1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1])
transcript2 = np.array([0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1])

## DEFINE DIMENSIONS + COORDINATES FOR ELEMENTS

page_width = 2000
page_height = 1324

title_y = 40
title_size = 40

panel_gap_y = title_size 

header_y = title_size + panel_gap_y
header_height = page_height / 6

exon_panel_x = page_width / 28
exon_panel_y = header_y + header_height + panel_gap_y
exon_panel_width = page_width - 2*exon_panel_x
exon_panel_height = page_height / 3

transcript_line_offset_x = exon_panel_x + exon_panel_width * 0.1
transcript_line_offset_y = exon_panel_y + exon_panel_height * 0.1
transcript_line_width = exon_panel_width * 0.8

transcript_text_x = exon_panel_x + (transcript_line_offset_x - exon_panel_x)/2
transcript_text_y = transcript_line_offset_y

coverage_panel_x = exon_panel_x
coverage_panel_y = exon_panel_y + exon_panel_height + panel_gap_y
coverage_panel_width = exon_panel_width
coverage_panel_height = exon_panel_height 
 


## CREATE DRAWER
d = drawer("test_output.png", width = page_width, height = page_height)

## CREATE HEADER
d.draw_text(text="SIRVsuite: coverage", y = title_y, x = page_width/2, h_align="center", font_size=title_size)

tab = {"Sample":"Sample 0001",
"Experiment":"Drosophila test experiment",
"Gene":"SIRV1",
"Mode":"complete coverage"}
d.draw_table(table_dict=tab,x=exon_panel_x,y=header_y,width=exon_panel_width/3,height=header_height)

## DRAW TRANSCRIPT PANEL

starts_in = np.array([d.return_differences(transcript1,mode="positive"),d.return_differences(transcript2,mode="positive")]).flatten()
ends_in = np.array([d.return_differences(transcript1,mode="negative"),d.return_differences(transcript2,mode="negative")]).flatten()

start, ends = get_continous_coverage_ends(starts_in,ends_in)
segment_lengths = [ends[i] - start[i] for i in range(len(start))]
total_segment_length = sum(segment_lengths)
intersegment_gap = 40
num_segments = len(start)
draw_length = transcript_line_width - (num_segments+1)*intersegment_gap

d.draw_rectangle(x = exon_panel_x,
                 y = exon_panel_y,
                 width = exon_panel_width,
                 height = exon_panel_height,
                 color_fill=(.5,.5,.5),
                 alpha = 0.2)

d.draw_line(x = transcript_line_offset_x,
            y = transcript_line_offset_y,
            width = transcript_line_width,
            line_width = 2, alpha=0.5)

segment_start = transcript_line_offset_x + intersegment_gap
for i in segment_lengths:
    d.draw_rectangle(x=segment_start, y=transcript_line_offset_y-20, width=i/total_segment_length*draw_length, height=40, alpha=0.8, color_fill=(60/255,140/255,80/255))
    segment_start += i/total_segment_length*draw_length + intersegment_gap

d.draw_text(text = t_name, 
            x = transcript_text_x, 
            y = transcript_text_y, 
            font_size = 28,
            v_align="center")


d.draw_text(text = "EXON DISTRIBUTION",x=exon_panel_x+120,y=exon_panel_y+exon_panel_height/2+140,font_size=28, v_align="center", h_align="center", rotate=270)

## DRAW COVERAGE PANEL
d.draw_rectangle(x = coverage_panel_x,
                 y = coverage_panel_y,
                 width = coverage_panel_width,
                 height = coverage_panel_height,
                 color_fill=(.5,.5,.5),
                 alpha = 0.2)

d.draw_text(text = "COVERAGE",x=exon_panel_x+60,y=coverage_panel_y+coverage_panel_height/2+80,font_size=28, v_align="center", h_align="center", rotate=270)

d.draw_signal(x = transcript_line_offset_x,
              y = coverage_panel_y,
              width = transcript_line_width,
              height = coverage_panel_height,
              signal = signal_normal,
              mode = "normal",
              alpha = 0.5,
              y_max = 51,
              color_fill = (.2,1,.6),
              line_width = 1,
              color_line = (0,.8,1),
              upside_down = False)

d.draw_signal(x = transcript_line_offset_x,
              y = coverage_panel_y,
              width = transcript_line_width,
              height = coverage_panel_height,
              signal = signal_segment,
              mode = "segment",
              alpha = 0.5,
              color_fill = (.1,.1,.1),
              line_width = 1,
              color_line = (0,.8,1),
              upside_down = False)

# draw y axis

d.draw_line(x=exon_panel_x+(transcript_line_offset_x-exon_panel_x)/2, y=coverage_panel_y+coverage_panel_height/2 - 10, width=coverage_panel_height/2-10, end_shape=("right","arrow"), rotate=270)
d.draw_line(x=exon_panel_x+(transcript_line_offset_x-exon_panel_x)/2, y=coverage_panel_y+coverage_panel_height/2 + 10, width=coverage_panel_height/2 - 10, end_shape=("right","arrow"), rotate=90)
d.draw_line(x=exon_panel_x+(transcript_line_offset_x-exon_panel_x)/2-10, y=coverage_panel_y+coverage_panel_height/2, width = 20)

# draw x axis

# draw result table
tab2 = {"CoD(+)":"0", "CoD(-)":"0.1231", "total_reads(+)": "4567", "total_reads(-)": "123"}
d.draw_rectangle(x=exon_panel_x+exon_panel_width*2/3-20,y=header_y-20,width=exon_panel_width/3+20,height=header_height+40, color_fill=(.5,.5,.5), alpha=.2)
d.draw_table(table_dict=tab2,x=exon_panel_x+exon_panel_width*2/3,y=header_y,width=exon_panel_width/3,height=header_height)


d.finish()