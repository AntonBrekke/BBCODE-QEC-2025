import src.BBcode_classes as BBcode_classes
import src.TileCode_classes as TileCode_classes
import src.decoder_classes as decoder_classes
from panqec.gui import GUI


"""
* Code written by Anton Brekke * 

This file runs all BBcodes and decoders made in the "src" folder, and visualizes it in PanQEC's GUI. 
To visualize the grid, a .json file is needed. The class "BB2DCode" makes use of _BBCode.json. 
"""

gui = GUI()

code_class_dict = {}
for name, cls in BBcode_classes.__dict__.items():
    if isinstance(cls, type) and issubclass(cls, BBcode_classes.BB2DCode) and name != 'BB2DCode':
        code_class_dict[name] = cls
        gui.add_code(cls, name)

code_class_dict_TileCode = {}
for name, cls in TileCode_classes.__dict__.items():
    if isinstance(cls, type) and issubclass(cls, TileCode_classes.TileCode) and name != 'TileCode':
        code_class_dict_TileCode[name] = cls
        gui.add_code(cls, name)

decoder_dict = {}
for name, cls in decoder_classes.__dict__.items():
    if isinstance(cls, type) and issubclass(cls, decoder_classes.BaseDecoder) and name != 'BaseDecoder':
        decoder_dict[name] = cls
        gui.add_decoder(cls, name)

gui.run(port=5000)