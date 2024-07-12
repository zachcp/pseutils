// Use this builder to navigate the creating of MolViewSpec files. Chain operations together as needed and invoke
//`get_state` or `save_state` to export the corresponding JSON needed to recreate that scene.

use crate::molviewspec::nodes::{DescriptionFormatT, KindT, Node, State};

pub fn create_builder() -> State {
    return State::new();
}

// struct Base {
//     root: PrivateAttr<Root>,
//     node: PrivateAttr<Node>,
// }

// impl Base {
//     fn new(root: Root, node: Node) -> Self {
//         Self {
//             root: PrivateAttr::new(root),
//             node: PrivateAttr::new(node),
//         }
//     }

//     fn add_child(&mut self, node: Node) {
//         if self.node.children.is_none() {
//             self.node.children = Some(Vec::new());
//         }
//         self.node.children.as_mut().unwrap().push(node);
//     }
// }

// // Manually position the camera.
// pub fn camera(
//     &mut self,
//     target: (f32, f32, f32),
//     position: (f32, f32, f32),
//     up: Option<(f32, f32, f32)>,
// ) -> &mut Self {
//     let params = CameraParams {
//         target,
//         position,
//         up: up.unwrap_or((0.0, 1.0, 0.0)),
//     };
//     let node = Node::new_with_params("camera", params);
//     self.add_child(node);
//     self
// }

// // Customize canvas properties such as background color.
// pub fn canvas(&mut self, background_color: Option<ColorT>) -> &mut Self {
//     let params = CanvasParams { background_color };
//     let node = Node::new_with_params("canvas", params);
//     self.add_child(node);
//     self
// }

// // Add a new structure to the builder by downloading structure data from a URL.
// pub fn download(&mut self, url: String) -> Download {
//     let params = DownloadParams { url };
//     let node = Node::new_with_params("download", params);
//     self.add_child(node.clone());
//     Download::new(node, self)
// }

// /// Experimental: Allows the definition of generic visuals such as spheres and lines.
// pub fn generic_visuals(&mut self) -> GenericVisuals {
//     let node = Node::new("generic_visuals");
//     self.add_child(node.clone());
//     GenericVisuals::new(node, self)
// }

// fn add_child(&mut self, node: Node) {
//     self.node.add_child(node);
// }
// }
//
// impl Download {
//     /// Builder step with operations needed after downloading structure data.
//     pub fn new(node: Node, root: Rc<RefCell<Root>>) -> Self {
//         Self {
//             base: Base::new(node, root),
//         }
//     }

//     /// Parse the content by specifying the file format.
//     pub fn parse(&mut self, format: ParseFormatT) -> Parse {
//         let params = make_params(ParseParams::new(format));
//         let node = Node::new("parse", params);
//         self.base.add_child(node.clone());
//         Parse::new(node, self.base.root.clone())
//     }
// }

// impl Parse {
//     /// Builder step with operations needed after parsing structure data.
//     pub fn new(node: Node, root: Rc<RefCell<Root>>) -> Self {
//         Self {
//             base: Base::new(node, root),
//         }
//     }

//     /// Create a structure for the deposited coordinates.
//     pub fn model_structure(
//         &mut self,
//         model_index: Option<i32>,
//         block_index: Option<i32>,
//         block_header: Option<String>,
//     ) -> Structure {
//         let params = make_params(StructureParams::new(
//             "model",
//             model_index,
//             block_index,
//             block_header,
//             None,
//             None,
//             None,
//         ));
//         let node = Node::new("structure", params);
//         self.base.add_child(node.clone());
//         Structure::new(node, self.base.root.clone())
//     }

//     /// Create an assembly structure.
//     pub fn assembly_structure(
//         &mut self,
//         assembly_id: Option<String>,
//         model_index: Option<i32>,
//         block_index: Option<i32>,
//         block_header: Option<String>,
//     ) -> Structure {
//         let params = make_params(StructureParams::new(
//             "assembly",
//             model_index,
//             block_index,
//             block_header,
//             assembly_id,
//             None,
//             None,
//         ));
//         let node = Node::new("structure", params);
//         self.base.add_child(node.clone());
//         Structure::new(node, self.base.root.clone())
//     }

//     /// Create symmetry structure for a given range of Miller indices.
//     pub fn symmetry_structure(
//         &mut self,
//         ijk_min: Option<(i32, i32, i32)>,
//         ijk_max: Option<(i32, i32, i32)>,
//         model_index: Option<i32>,
//         block_index: Option<i32>,
//         block_header: Option<String>,
//     ) -> Structure {
//         let params = make_params(StructureParams::new(
//             "symmetry",
//             model_index,
//             block_index,
//             block_header,
//             None,
//             ijk_min.or(Some((-1, -1, -1))),
//             ijk_max.or(Some((1, 1, 1))),
//         ));
//         let node = Node::new("structure", params);
//         self.base.add_child(node.clone());
//         Structure::new(node, self.base.root.clone())
//     }

//     /// Create structure of symmetry mates.
//     pub fn symmetry_mates_structure(
//         &mut self,
//         radius: Option<f32>,
//         model_index: Option<i32>,
//         block_index: Option<i32>,
//         block_header: Option<String>,
//     ) -> Structure {
//         let params = make_params(StructureParams::new(
//             "symmetry_mates",
//             model_index,
//             block_index,
//             block_header,
//             None,
//             None,
//             None,
//         ));
//         let node = Node::new("structure", params);
//         self.base.add_child(node.clone());
//         Structure::new(node, self.base.root.clone())
//     }
// }

// pub struct Structure {
//     base: Base,
// }

// impl Structure {
//     /// Builder step with operations needed after defining the structure to work with.
//     pub fn new(node: Node, root: Rc<RefCell<Root>>) -> Self {
//         Self {
//             base: Base::new(node, root),
//         }
//     }

//     /// Define a new component/selection for the given structure.
//     pub fn component(&mut self, selector: ComponentSelector) -> Component {
//         let params = make_params(ColorInlineParams::new(selector));
//         let node = Node::new("component", params);
//         self.base.add_child(node.clone());
//         Component::new(node, self.base.root.clone())
//     }

//     /// Define a new component/selection for the given structure by fetching additional data from a resource.
//     pub fn component_from_uri(
//         &mut self,
//         uri: String,
//         format: SchemaFormatT,
//         category_name: Option<String>,
//         field_name: Option<String>,
//         block_header: Option<String>,
//         block_index: Option<i32>,
//         schema: SchemaT,
//         field_values: Option<Vec<String>>,
//     ) -> Component {
//         let params = make_params(ComponentFromUriParams::new(
//             uri,
//             format,
//             category_name,
//             field_name,
//             block_header,
//             block_index,
//             schema,
//             field_values,
//         ));
//         let node = Node::new("component_from_uri", params);
//         self.base.add_child(node.clone());
//         Component::new(node, self.base.root.clone())
//     }

//     /// Define a new component/selection for the given structure by using categories from the source file.
//     pub fn component_from_source(
//         &mut self,
//         category_name: String,
//         field_name: Option<String>,
//         block_header: Option<String>,
//         block_index: Option<i32>,
//         schema: SchemaT,
//         field_values: Option<Vec<String>>,
//     ) -> Component {
//         let params = make_params(ComponentFromSourceParams::new(
//             category_name,
//             field_name,
//             block_header,
//             block_index,
//             schema,
//             field_values,
//         ));
//         let node = Node::new("component_from_source", params);
//         self.base.add_child(node.clone());
//         Component::new(node, self.base.root.clone())
//     }

//     /// Define a new label for the given structure by fetching additional data from a resource.
//     pub fn label_from_uri(
//         &mut self,
//         uri: String,
//         format: SchemaFormatT,
//         category_name: Option<String>,
//         field_name: Option<String>,
//         block_header: Option<String>,
//         block_index: Option<i32>,
//         schema: SchemaT,
//     ) -> &mut Self {
//         let params = make_params(LabelFromUriParams::new(
//             uri,
//             format,
//             category_name,
//             field_name,
//             block_header,
//             block_index,
//             schema,
//         ));
//         let node = Node::new("label_from_uri", params);
//         self.base.add_child(node);
//         self
//     }

//     /// Define a new label for the given structure by fetching additional data from the source file.
//     pub fn label_from_source(
//         &mut self,
//         category_name: String,
//         field_name: Option<String>,
//         block_header: Option<String>,
//         block_index: Option<i32>,
//         schema: SchemaT,
//     ) -> &mut Self {
//         let params = make_params(LabelFromSourceParams::new(
//             category_name,
//             field_name,
//             block_header,
//             block_index,
//             schema,
//         ));
//         let node = Node::new("label_from_source", params);
//         self.base.add_child(node);
//         self
//     }

//     /// Define a new tooltip for the given structure by fetching additional data from a resource.
//     pub fn tooltip_from_uri(
//         &mut self,
//         uri: String,
//         format: SchemaFormatT,
//         category_name: Option<String>,
//         field_name: Option<String>,
//         block_header: Option<String>,
//         block_index: Option<i32>,
//         schema: SchemaT,
//     ) -> &mut Self {
//         let params = make_params(TooltipFromUriParams::new(
//             uri,
//             format,
//             category_name,
//             field_name,
//             block_header,
//             block_index,
//             schema,
//         ));
//         let node = Node::new("tooltip_from_uri", params);
//         self.base.add_child(node);
//         self
//     }

//     /// Define a new tooltip for the given structure by fetching additional data from the source file.
//     pub fn tooltip_from_source(
//         &mut self,
//         category_name: String,
//         field_name: Option<String>,
//         block_header: Option<String>,
//         block_index: Option<i32>,
//         schema: SchemaT,
//     ) -> &mut Self {
//         let params = make_params(TooltipFromSourceParams::new(
//             category_name,
//             field_name,
//             block_header,
//             block_index,
//             schema,
//         ));
//         let node = Node::new("tooltip_from_source", params);
//         self.base.add_child(node);
//         self
//     }

//     /// Transform a structure by applying a rotation matrix and/or translation vector.
//     pub fn transform(
//         &mut self,
//         rotation: Option<[f32; 9]>,
//         translation: Option<[f32; 3]>,
//     ) -> Result<&mut Self, String> {
//         if let Some(rot) = rotation {
//             if !Self::is_rotation_matrix(&rot) {
//                 return Err("Parameter `rotation` must be a rotation matrix".to_string());
//             }
//         }

//         let params = make_params(TransformParams::new(rotation, translation));
//         let node = Node::new("transform", params);
//         self.base.add_child(node);
//         Ok(self)
//     }

//     fn is_rotation_matrix(t: &[f32; 9], eps: f32 = 0.005) -> bool {
//         let [a00, a01, a02, a10, a11, a12, a20, a21, a22] = *t;

//         let det3x3 = (a00 * (a11 * a22 - a12 * a21) - a01 * (a10 * a22 - a12 * a20) + a02 * (a10 * a21 - a11 * a20)).abs();
//         (det3x3 - 1.0).abs() < eps
//     }
// }

// pub struct Component {
//     base: Base,
// }

// impl Component {
//     /// Builder step with operations relevant for a particular component.
//     pub fn new(node: Node, root: Rc<RefCell<Root>>) -> Self {
//         Self {
//             base: Base::new(node, root),
//         }
//     }

//     /// Add a representation for this component.
//     pub fn representation(&mut self, rep_type: RepresentationTypeT) -> Representation {
//         let params = make_params(RepresentationParams::new(rep_type));
//         let node = Node::new("representation", params);
//         self.base.add_child(node.clone());
//         Representation::new(node, self.base.root.clone())
//     }

//     /// Add a text label to a component.
//     pub fn label(&mut self, text: String) -> &mut Self {
//         let params = make_params(LabelInlineParams::new(text));
//         let node = Node::new("label", params);
//         self.base.add_child(node);
//         self
//     }

//     /// Add a tooltip that shows additional information of a component when hovering over it.
//     pub fn tooltip(&mut self, text: String) -> &mut Self {
//         let params = make_params(TooltipInlineParams::new(text));
//         let node = Node::new("tooltip", params);
//         self.base.add_child(node);
//         self
//     }

//     /// Focus on this structure or component.
//     pub fn focus(
//         &mut self,
//         direction: Option<(f32, f32, f32)>,
//         up: Option<(f32, f32, f32)>,
//     ) -> &mut Self {
//         let params = make_params(FocusInlineParams::new(direction, up));
//         let node = Node::new("focus", params);
//         self.base.add_child(node);
//         self
//     }
// }

// struct Base {
//     node: Node,
//     root: Rc<RefCell<Root>>,
// }

// impl Base {
//     fn new(node: Node, root: Rc<RefCell<Root>>) -> Self {
//         Self { node, root }
//     }

//     fn add_child(&mut self, node: Node) {
//         self.node.add_child(node);
//     }
// }

// impl Representation {
//     /// Builder step with operations relating to particular representations.
//     pub fn new(node: Node, root: Rc<RefCell<Root>>) -> Self {
//         Self {
//             base: Base::new(node, root),
//         }
//     }

//     /// Use a custom category from the source file to define colors of this representation.
//     pub fn color_from_source(
//         &mut self,
//         schema: SchemaT,
//         category_name: String,
//         field_name: Option<String>,
//         block_header: Option<String>,
//         block_index: Option<i32>,
//     ) -> &mut Self {
//         let params = make_params(ColorFromSourceParams::new(
//             schema,
//             category_name,
//             field_name,
//             block_header,
//             block_index,
//         ));
//         let node = Node::new("color_from_source", params);
//         self.base.add_child(node);
//         self
//     }

//     /// Use another resource to define colors of this representation.
//     pub fn color_from_uri(
//         &mut self,
//         schema: SchemaT,
//         uri: String,
//         format: String,
//         category_name: Option<String>,
//         field_name: Option<String>,
//         block_header: Option<String>,
//         block_index: Option<i32>,
//     ) -> &mut Self {
//         let params = make_params(ColorFromUriParams::new(
//             schema,
//             uri,
//             format,
//             category_name,
//             field_name,
//             block_header,
//             block_index,
//         ));
//         let node = Node::new("color_from_uri", params);
//         self.base.add_child(node);
//         self
//     }

//     /// Customize the color of this representation.
//     pub fn color(&mut self, color: ColorT, selector: ComponentSelector) -> &mut Self {
//         let params = make_params(ColorInlineParams::new(color, selector));
//         let node = Node::new("color", params);
//         self.base.add_child(node);
//         self
//     }
// }

// pub struct GenericVisuals {
//     base: Base,
// }

// impl GenericVisuals {
//     /// Experimental builder for custom, primitive visuals.
//     pub fn new(node: Node, root: Rc<RefCell<Root>>) -> Self {
//         Self {
//             base: Base::new(node, root),
//         }
//     }

//     /// Draw a sphere.
//     pub fn sphere(
//         &mut self,
//         position: (f32, f32, f32),
//         radius: f32,
//         color: ColorT,
//         label: Option<String>,
//         tooltip: Option<String>,
//     ) -> &mut Self {
//         let params = make_params(SphereParams::new(position, radius, color, label, tooltip));
//         let node = Node::new("sphere", params);
//         self.base.add_child(node);
//         self
//     }

//     /// Draw a line.
//     pub fn line(
//         &mut self,
//         position1: (f32, f32, f32),
//         position2: (f32, f32, f32),
//         radius: f32,
//         color: ColorT,
//         label: Option<String>,
//         tooltip: Option<String>,
//     ) -> &mut Self {
//         let params = make_params(LineParams::new(
//             position1, position2, radius, color, label, tooltip,
//         ));
//         let node = Node::new("line", params);
//         self.base.add_child(node);
//         self
//     }
// }

// fn make_params<T>(params: T) -> T {
//     params
// }
