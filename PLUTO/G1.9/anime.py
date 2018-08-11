# Recorded script from Mayavi2
from numpy import array
try:
    engine = mayavi.engine
except NameError:
    from mayavi.api import Engine
    engine = Engine()
    engine.start()
if len(engine.scenes) == 0:
    engine.new_scene()
# ------------------------------------------- 
scene = engine.scenes[0]
scene.scene.camera.position = [311.77701153000146, 311.77701153000146, 311.77701153000146]
scene.scene.camera.focal_point = [64.5, 64.5, 64.5]
scene.scene.camera.view_angle = 30.0
scene.scene.camera.view_up = [0.0, 0.0, 1.0]
scene.scene.camera.clipping_range = [203.4193931866478, 712.40317819590632]
scene.scene.camera.compute_view_plane_normal()
scene.scene.render()
image_plane_widget1 = engine.scenes[0].children[1].children[0].children[0]
image_plane_widget1.ipw.origin = array([ 0.5,  1. ,  0.5])
image_plane_widget1.ipw.slice_index = 0
image_plane_widget1.ipw.slice_position = 1.0
image_plane_widget1.ipw.point1 = array([   0.5,    1. ,  128.5])
image_plane_widget1.ipw.point2 = array([ 128.5,    1. ,    0.5])
image_plane_widget1.ipw.origin = array([ 0.5,  1. ,  0.5])
image_plane_widget1.ipw.point1 = array([   0.5,    1. ,  128.5])
image_plane_widget1.ipw.point2 = array([ 128.5,    1. ,    0.5])
image_plane_widget2 = engine.scenes[0].children[2].children[0].children[0]
image_plane_widget2.ipw.origin = array([ 0.5,  0.5,  1. ])
image_plane_widget2.ipw.slice_index = 0
image_plane_widget2.ipw.slice_position = 1.0
image_plane_widget2.ipw.point1 = array([ 128.5,    0.5,    1. ])
image_plane_widget2.ipw.point2 = array([   0.5,  128.5,    1. ])
image_plane_widget2.ipw.origin = array([ 0.5,  0.5,  1. ])
image_plane_widget2.ipw.point1 = array([ 128.5,    0.5,    1. ])
image_plane_widget2.ipw.point2 = array([   0.5,  128.5,    1. ])
image_plane_widget = engine.scenes[0].children[0].children[0].children[0]
image_plane_widget.ipw.origin = array([ 1. ,  0.5,  0.5])
image_plane_widget.ipw.slice_index = 0
image_plane_widget.ipw.slice_position = 1.0
image_plane_widget.ipw.point1 = array([   1. ,  128.5,    0.5])
image_plane_widget.ipw.point2 = array([   1. ,    0.5,  128.5])
image_plane_widget.ipw.origin = array([ 1. ,  0.5,  0.5])
image_plane_widget.ipw.point1 = array([   1. ,  128.5,    0.5])
image_plane_widget.ipw.point2 = array([   1. ,    0.5,  128.5])
scene.scene.show_axes = True
scene.scene.camera.position = [486.97762996578587, 32.807552030816318, 127.31721271939678]
scene.scene.camera.focal_point = [64.5, 64.5, 64.5]
scene.scene.camera.view_angle = 30.0
scene.scene.camera.view_up = [-0.16242997732460931, -0.57282261680926339, 0.80342439105252128]
scene.scene.camera.clipping_range = [270.27987283787792, 628.23960958969758]
scene.scene.camera.compute_view_plane_normal()
scene.scene.render()
from mayavi.core.module_manager import ModuleManager
module_manager3 = ModuleManager()
engine.add_module(module_manager3, obj=None)
module_manager3.scalar_lut_manager.reverse_lut = True
module_manager3.scalar_lut_manager.reverse_lut = False
module_manager3.scalar_lut_manager.reverse_lut = True
module_manager3.scalar_lut_manager.reverse_lut = False
array_source2 = engine.scenes[0].children[2]
array_source2.children[1:2] = []

module_manager4 = ModuleManager()
engine.add_module(module_manager4, obj=None)
array_source = engine.scenes[0].children[0]
array_source.children[1:2] = []
from mayavi.modules.axes import Axes
axes = Axes()
module_manager = engine.scenes[0].children[0].children[0]
engine.add_filter(axes, module_manager)
scene.scene.camera.position = [437.19923996106922, -121.22299589742653, 164.6998334377177]
scene.scene.camera.focal_point = [64.5, 64.5, 64.5]
scene.scene.camera.view_angle = 30.0
scene.scene.camera.view_up = [-0.4727384608678496, -0.60095494152678108, 0.64449321631095513]
scene.scene.camera.clipping_range = [228.1627604400266, 681.25637670357492]
scene.scene.camera.compute_view_plane_normal()
scene.scene.render()
from mayavi.modules.surface import Surface
surface = Surface()
engine.add_filter(surface, module_manager)
scene.scene.camera.position = [457.54060862037477, 151.6004944233332, -81.685995020579185]
scene.scene.camera.focal_point = [64.5, 64.5, 64.5]
scene.scene.camera.view_angle = 30.0
scene.scene.camera.view_up = [0.37956817215096489, -0.70260197392637047, 0.60189572927936885]
scene.scene.camera.clipping_range = [237.76612866853119, 669.16771468729144]
scene.scene.camera.compute_view_plane_normal()
scene.scene.render()
scene.scene.camera.position = [384.2249016078232, 347.25880978783448, 100.01343510357127]
scene.scene.camera.focal_point = [64.5, 64.5, 64.5]
scene.scene.camera.view_angle = 30.0
scene.scene.camera.view_up = [0.48118020237428633, -0.62170730472395053, 0.61801750792022203]
scene.scene.camera.clipping_range = [234.29587407191352, 673.53605024735543]
scene.scene.camera.compute_view_plane_normal()
scene.scene.render()
scene.scene.camera.position = [467.53271647812761, 86.862090969786067, 207.68633883628965]
scene.scene.camera.focal_point = [64.5, 64.5, 64.5]
scene.scene.camera.view_angle = 30.0
scene.scene.camera.view_up = [-0.13964044259461819, -0.84004030228192661, 0.5242450165085677]
scene.scene.camera.clipping_range = [254.93770732290054, 647.55223502437934]
scene.scene.camera.compute_view_plane_normal()
scene.scene.render()
scene.scene.camera.position = [488.8868933754232, 12.448309188877186, 89.482944351718913]
scene.scene.camera.focal_point = [64.5, 64.5, 64.5]
scene.scene.camera.view_angle = 30.0
scene.scene.camera.view_up = [-0.13400023770542915, -0.84074597935105932, 0.52458567889327357]
scene.scene.camera.clipping_range = [274.90857320997924, 622.41302947305701]
scene.scene.camera.compute_view_plane_normal()
scene.scene.render()
# ------------------------------------------- 
from mayavi.tools.show import show
show()
