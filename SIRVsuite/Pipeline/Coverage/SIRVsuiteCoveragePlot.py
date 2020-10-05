class SIRVsuiteCoveragePlot():
    def __init__(self):
        self.page_width = 2000
        self.page_height = 1324

        self.exon_panel_x = page_width / 28
        self.exon_panel_y = page_height / 4
        self.exon_panel_width = page_width - 2*exon_panel_x
        self.exon_panel_height = page_height / 3

        transcript_line_offset_x = exon_panel_x + exon_panel_width * 0.1
        transcript_line_offset_y = exon_panel_y + exon_panel_height * 0.1
        transcript_line_width = exon_panel_width * 0.8

        transcript_text_x = exon_panel_x + (transcript_line_offset_x - exon_panel_x)/2
        transcript_text_y = transcript_line_offset_y

        panel_gap_y = exon_panel_y / 8

        coverage_panel_x = exon_panel_x
        coverage_panel_y = exon_panel_y + exon_panel_height + panel_gap_y
        coverage_panel_width = exon_panel_width
        coverage_panel_height = exon_panel_height 
    """
    def draw_header(self, info_dict):

    def draw_statistics_table():
        

    def draw_transcript_panel(self, exon_coords):

    def draw_coverage_panel(self, real_coverage):

    def get_continous_coverage_ends(starts, ends):
    
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
        """
