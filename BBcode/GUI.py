import source.BBcode_classes as BBcode_classes
import source.decoder_classes as decoder_classes
from panqec.gui import GUI


gui = GUI()

code_class_dict = {}
for name, cls in BBcode_classes.__dict__.items():
    if isinstance(cls, type) and issubclass(cls, BBcode_classes.AntonBB2DCode) and name != 'AntonBB2DCode':
        code_class_dict[name] = cls
        gui.add_code(cls, name)

decoder_dict = {}
for name, cls in decoder_classes.__dict__.items():
    if isinstance(cls, type) and issubclass(cls, decoder_classes.BaseDecoder) and name != 'BaseDecoder':
        decoder_dict[name] = cls
        gui.add_decoder(cls, name)

gui.run(port=5000)