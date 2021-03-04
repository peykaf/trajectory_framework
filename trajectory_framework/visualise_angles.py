import numpy
from traits.api import HasTraits, Range, Instance, on_trait_change
from traitsui.api import View, Item, Group
from tvtk.pyface.scene_editor import SceneEditor
from mayavi.tools.mlab_scene_model import MlabSceneModel
from mayavi.core.api import PipelineBase
from mayavi.core.ui.mayavi_scene import MayaviScene

from mayavi.mlab import points3d, quiver3d
from scipy.spatial import ConvexHull, Delaunay

class TrajectoryVisualization(HasTraits):
    # the layout of the dialog created
    scene = Instance(MlabSceneModel, ())
    plot = Instance(PipelineBase)
    control_point = Range(0, 103, 0)
    # The layout of the dialog created
    view = View(Item('scene', editor=SceneEditor(scene_class=MayaviScene),
                     height=500, width=600, show_label=False),
                Group(
                        '_', 'angle',
                     ),
                resizable=True,
                )

    def __init__(self, couch_points):
        # Do not forget to call the parent's __init__
        HasTraits.__init__(self)

        self.couch_points = couch_points
        self.setup_animation()

    @on_trait_change('angle')
    def update_plot(self):
        self.change_precalculated(update=True)

    def change_precalculated(self, update=False):
        print "Theta: %f, phi: %f" % (self.control_points[self.control_point]["gantry_angle"] * 180 / numpy.pi, self.control_points[self.control_point]["couch_angle"] * 180 / numpy.pi)
        x_dirs = self.control_points[self.angle]["x_dirs"]
        y_dirs = self.control_points[self.angle]["y_dirs"]
        z_dirs = self.control_points[self.angle]["z_dirs"]

        x_pos = self.control_points[self.angle]["x_pos"]
        y_pos = self.control_points[self.angle]["y_pos"]
        z_pos = self.control_points[self.angle]["z_pos"]

        x_couch = self.control_points[self.angle]["x_couch"]
        y_couch = self.control_points[self.angle]["y_couch"]
        z_couch = self.control_points[self.angle]["z_couch"]

        if update:
            self.plot.mlab_source.set(x=x_pos, y=y_pos, z=z_pos, u=x_dirs, v=y_dirs, w=z_dirs, scale_factor=100)
            
            new_x_points = [bla[0] + (x_couch - self.first_couch[0]) for bla in self.couch_triangles]
            new_y_points = [bla[1] + (y_couch - self.first_couch[1]) for bla in self.couch_triangles]
            new_z_points = [bla[2] + (z_couch - self.first_couch[2]) for bla in self.couch_triangles]

            self.couch.mlab_source.set(x=new_x_points, y=new_y_points, z=new_z_points)

        else:
            self.plot = self.scene.mlab.quiver3d(x_pos, y_pos, z_pos, x_dirs, y_dirs, z_dirs, scalars=scores, scale_factor=100, colormap="hot")

    def visualise_couch(self):
        first_x = self.control_points[0]["x_couch"]
        first_y = self.control_points[0]["y_couch"]
        first_z = self.control_points[0]["z_couch"]
        self.first_couch = [first_x, first_y, first_z]

        couch_pos = []
        for point in self.couch_points:
            couch_pos.append([pos + self.first_couch[i] for i, pos in enumerate(point)])

        couch_convex_hull = ConvexHull(couch_pos)
        self.couch_triangles = couch_convex_hull.points[couch_convex_hull.simplices]
        x_points = []
        y_points = []
        z_points = []
        for triangle in self.couch_triangles:
            x_points += [bla[0] for bla in triangle]
            y_points += [bla[1] for bla in triangle]
            z_points += [bla[2] for bla in triangle]

        triangles = [[3*n, 3*n+1, 3*n+2] for n in range(len(x_points) / 3)]

        self.couch = self.scene.mlab.triangular_mesh(x_points, y_points, z_points, triangles, color=(1, 0.9, 0.8), opacity=1.0)

    def setup_animation(self):
        self.visualise_couch()
        self.change_precalculated(update=False)


