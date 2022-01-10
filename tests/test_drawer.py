from _pytest.fixtures import fixture
from numpy.core.fromnumeric import round_
import pytest

from SIRVsuite.Pipeline.Coverage.cairoDrawer import CairoDrawer as cd
import shutil
import os

@pytest.mark.parametrize("export_types", ["svg", "png", "ps", "pdf"])

def test_drawer(export_types):
    drawer = cd(out_path="./tests/test_data/temp/test_graphics."+export_types)
    
    # Line tests
    drawer.draw_line(x=10, y=10, width=100, end_shape=("both", "line"))
    drawer.draw_line(x=10, y=30, width=100, end_shape=("both", "arrow"))
    drawer.draw_line(x=10, y=50, width=100)
    
    # Rectangle tests
    drawer.draw_rectangle(x=70, y=70, width=200, height=200, round_aspect=0.8)
    
    # Text tests
    drawer.draw_text(text="TEST", x=100, y=400, h_align="left")
    drawer.draw_text(text="TEST", x=100, y=450, h_align="right")
    drawer.draw_text(text="TEST", x=100, y=500, h_align="center")
    
    # Signal tests
    drawer.draw_signal(signal=[1,2,3,4,5,4,3,2,1], x=100, y=700, width=200, height=100, y_max=0)
    drawer.draw_signal(signal=[1,2,3,4,5,4,3,2,1], x=100, y=900, width=200, height=100)
    drawer.draw_signal(signal=[0,1,2,3,4,5,6,7,8,9,10,9,8,7,6,5,4,3,2,1,0], x=100, y=1100, width=200, height=100, mode="segment")
    
    drawer.finish()

    if os.path.exists("./tests/test_data/temp"):
       shutil.rmtree("./tests/test_data/temp")
