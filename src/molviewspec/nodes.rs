use serde::{Deserialize, Serialize};
use std::collections::HashMap;

#[derive(Serialize, Deserialize, Debug)]
#[serde(rename_all = "snake_case")]
enum KindT {
    Root,
    Camera,
    Canvas,
    Color,
    ColorFromSource,
    ColorFromUri,
    Component,
    ComponentFromSource,
    ComponentFromUri,
    Download,
    Focus,
    GenericVisuals,
    Label,
    LabelFromSource,
    LabelFromUri,
    Line,
    Parse,
    Representation,
    Sphere,
    Structure,
    Tooltip,
    TooltipFromSource,
    TooltipFromUri,
    Transform,
}

#[derive(Serialize, Deserialize, Debug)]
struct Node {
    kind: KindT,
    #[serde(skip_serializing_if = "Option::is_none")]
    params: Option<HashMap<String, serde_json::Value>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    children: Option<Vec<Node>>,
}

#[derive(Serialize, Deserialize, Debug)]
#[serde(rename_all = "snake_case")]
enum DescriptionFormatT {
    Markdown,
    Plaintext,
}

#[derive(Serialize, Deserialize, Debug)]
struct Metadata {
    version: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    title: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    description: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    description_format: Option<DescriptionFormatT>,
    timestamp: String,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct State {
    root: Node,
    metadata: Metadata,
}

#[derive(Serialize, Deserialize, Debug)]
struct DownloadParams {
    url: String,
}

#[derive(Serialize, Deserialize, Debug)]
#[serde(rename_all = "lowercase")]
enum ParseFormatT {
    Mmcif,
    Bcif,
    Pdb,
}

#[derive(Serialize, Deserialize, Debug)]
struct ParseParams {
    format: ParseFormatT,
}

#[derive(Serialize, Deserialize, Debug)]
#[serde(rename_all = "snake_case")]
enum StructureTypeT {
    Model,
    Assembly,
    Symmetry,
    SymmetryMates,
}

#[derive(Serialize, Deserialize, Debug)]
struct StructureParams {
    #[serde(rename = "type")]
    structure_type: StructureTypeT,
    #[serde(skip_serializing_if = "Option::is_none")]
    assembly_id: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    assembly_index: Option<i32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    model_index: Option<i32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    block_index: Option<i32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    block_header: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    radius: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    ijk_min: Option<(i32, i32, i32)>,
    #[serde(skip_serializing_if = "Option::is_none")]
    ijk_max: Option<(i32, i32, i32)>,
}

#[derive(Serialize, Deserialize, Debug)]
#[serde(rename_all = "snake_case")]
enum ComponentSelectorT {
    All,
    Polymer,
    Protein,
    Nucleic,
    Branched,
    Ligand,
    Ion,
    Water,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct ComponentExpression {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub label_entity_id: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub label_asym_id: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub auth_asym_id: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub label_seq_id: Option<i32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub auth_seq_id: Option<i32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub pdbx_pdb_ins_code: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub beg_label_seq_id: Option<i32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub end_label_seq_id: Option<i32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub beg_auth_seq_id: Option<i32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub end_auth_seq_id: Option<i32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub residue_index: Option<i32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub label_atom_id: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub auth_atom_id: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub type_symbol: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub atom_id: Option<i32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub atom_index: Option<i32>,
}

#[derive(Serialize, Deserialize, Debug)]
#[serde(rename_all = "snake_case")]
enum RepresentationTypeT {
    BallAndStick,
    Cartoon,
    Surface,
}

#[derive(Serialize, Deserialize, Debug)]
#[serde(rename_all = "lowercase")]
enum ColorNamesT {
    Aliceblue,
    Antiquewhite,
    Aqua,
    Aquamarine,
    Azure,
    Beige,
    Bisque,
    Black,
    Blanchedalmond,
    Blue,
    Blueviolet,
    Brown,
    Burlywood,
    Cadetblue,
    Chartreuse,
    Chocolate,
    Coral,
    Cornflowerblue,
    Cornsilk,
    Crimson,
    Cyan,
    Darkblue,
    Darkcyan,
    Darkgoldenrod,
    Darkgray,
    Darkgreen,
    Darkgrey,
    Darkkhaki,
    Darkmagenta,
    Darkolivegreen,
    Darkorange,
    Darkorchid,
    Darkred,
    Darksalmon,
    Darkseagreen,
    Darkslateblue,
    Darkslategray,
    Darkslategrey,
    Darkturquoise,
    Darkviolet,
    Deeppink,
    Deepskyblue,
    Dimgray,
    Dimgrey,
    Dodgerblue,
    Firebrick,
    Floralwhite,
    Forestgreen,
    Fuchsia,
    Gainsboro,
    Ghostwhite,
    Gold,
    Goldenrod,
    Gray,
    Green,
    Greenyellow,
    Grey,
    Honeydew,
    Hotpink,
    Indianred,
    Indigo,
    Ivory,
    Khaki,
    Lavender,
    Lavenderblush,
    Lawngreen,
    Lemonchiffon,
    Lightblue,
    Lightcoral,
    Lightcyan,
    Lightgoldenrodyellow,
    Lightgray,
    Lightgreen,
    Lightgrey,
    Lightpink,
    Lightsalmon,
    Lightseagreen,
    Lightskyblue,
    Lightslategray,
    Lightslategrey,
    Lightsteelblue,
    Lightyellow,
    Lime,
    Limegreen,
    Linen,
    Magenta,
    Maroon,
    Mediumaquamarine,
    Mediumblue,
    Mediumorchid,
    Mediumpurple,
    Mediumseagreen,
    Mediumslateblue,
    Mediumspringgreen,
    Mediumturquoise,
    Mediumvioletred,
    Midnightblue,
    Mintcream,
    Mistyrose,
    Moccasin,
    Navajowhite,
    Navy,
    Oldlace,
    Olive,
    Olivedrab,
    Orange,
    Orangered,
    Orchid,
    Palegoldenrod,
    Palegreen,
    Paleturquoise,
    Palevioletred,
    Papayawhip,
    Peachpuff,
    Peru,
    Pink,
    Plum,
    Powderblue,
    Purple,
    Red,
    Rosybrown,
    Royalblue,
    Saddlebrown,
    Salmon,
    Sandybrown,
    Seagreen,
    Seashell,
    Sienna,
    Silver,
    Skyblue,
    Slateblue,
    Slategray,
    Slategrey,
    Snow,
    Springgreen,
    Steelblue,
    Tan,
    Teal,
    Thistle,
    Tomato,
    Turquoise,
    Violet,
    Wheat,
    White,
    Whitesmoke,
    Yellow,
    Yellowgreen,
}

#[derive(Serialize, Deserialize, Debug)]
#[serde(untagged)]
enum ColorT {
    Named(ColorNamesT),
    Hex(String),
}

#[derive(Serialize, Deserialize, Debug)]
struct RepresentationParams {
    #[serde(rename = "type")]
    representation_type: RepresentationTypeT,
}

#[derive(Serialize, Deserialize, Debug)]
#[serde(rename_all = "snake_case")]
enum SchemaT {
    WholeStructure,
    Entity,
    Chain,
    AuthChain,
    Residue,
    AuthResidue,
    ResidueRange,
    AuthResidueRange,
    Atom,
    AuthAtom,
    AllAtomic,
}

#[derive(Serialize, Deserialize, Debug)]
#[serde(rename_all = "lowercase")]
enum SchemaFormatT {
    Cif,
    Bcif,
    Json,
}

#[derive(Serialize, Deserialize, Debug)]
struct DataFromUriParams {
    uri: String,
    format: SchemaFormatT,
    #[serde(skip_serializing_if = "Option::is_none")]
    category_name: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    field_name: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    block_header: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    block_index: Option<i32>,
    #[serde(rename = "schema")]
    schema_: SchemaT,
}

#[derive(Serialize, Deserialize, Debug)]
struct DataFromSourceParams {
    category_name: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    field_name: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    block_header: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    block_index: Option<i32>,
    #[serde(rename = "schema")]
    schema_: SchemaT,
}

#[derive(Serialize, Deserialize, Debug)]
struct ComponentInlineParams {
    selector: ComponentSelector,
}

#[derive(Serialize, Deserialize, Debug)]
#[serde(untagged)]
enum ComponentSelector {
    Selector(ComponentSelectorT),
    Expression(ComponentExpression),
    ExpressionList(Vec<ComponentExpression>),
}

#[derive(Serialize, Deserialize, Debug)]
struct ComponentFromUriParams {
    #[serde(flatten)]
    base: DataFromUriParams,
    #[serde(skip_serializing_if = "Option::is_none")]
    field_values: Option<Vec<String>>,
}

#[derive(Serialize, Deserialize, Debug)]
struct ComponentFromSourceParams {
    #[serde(flatten)]
    base: DataFromSourceParams,
    #[serde(skip_serializing_if = "Option::is_none")]
    field_values: Option<Vec<String>>,
}

#[derive(Serialize, Deserialize, Debug)]
struct ColorInlineParams {
    #[serde(flatten)]
    base: ComponentInlineParams,
    color: ColorT,
}

#[derive(Serialize, Deserialize, Debug)]
struct ColorFromUriParams {
    #[serde(flatten)]
    base: DataFromUriParams,
}

#[derive(Serialize, Deserialize, Debug)]
struct ColorFromSourceParams {
    #[serde(flatten)]
    base: DataFromSourceParams,
}

#[derive(Serialize, Deserialize, Debug)]
struct LabelInlineParams {
    text: String,
}

#[derive(Serialize, Deserialize, Debug)]
struct LabelFromUriParams {
    #[serde(flatten)]
    base: DataFromUriParams,
}

#[derive(Serialize, Deserialize, Debug)]
struct LabelFromSourceParams {
    #[serde(flatten)]
    base: DataFromSourceParams,
}

#[derive(Serialize, Deserialize, Debug)]
struct TooltipInlineParams {
    text: String,
}

#[derive(Serialize, Deserialize, Debug)]
struct TooltipFromUriParams {
    #[serde(flatten)]
    base: DataFromUriParams,
}

#[derive(Serialize, Deserialize, Debug)]
struct TooltipFromSourceParams {
    #[serde(flatten)]
    base: DataFromSourceParams,
}

#[derive(Serialize, Deserialize, Debug)]
struct FocusInlineParams {
    #[serde(skip_serializing_if = "Option::is_none")]
    direction: Option<(f32, f32, f32)>,
    #[serde(skip_serializing_if = "Option::is_none")]
    up: Option<(f32, f32, f32)>,
}

#[derive(Serialize, Deserialize, Debug)]
struct TransformParams {
    #[serde(skip_serializing_if = "Option::is_none")]
    rotation: Option<Vec<f32>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    translation: Option<(f32, f32, f32)>,
}

#[derive(Serialize, Deserialize, Debug)]
struct CameraParams {
    target: (f32, f32, f32),
    position: (f32, f32, f32),
    up: (f32, f32, f32),
}

#[derive(Serialize, Deserialize, Debug)]
struct CanvasParams {
    background_color: ColorT,
}

#[derive(Serialize, Deserialize, Debug)]
struct SphereParams {
    position: (f32, f32, f32),
    radius: f32,
    color: ColorT,
    #[serde(skip_serializing_if = "Option::is_none")]
    label: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    tooltip: Option<String>,
}

#[derive(Serialize, Deserialize, Debug)]
struct LineParams {
    position1: (f32, f32, f32),
    position2: (f32, f32, f32),
    radius: f32,
    color: ColorT,
    #[serde(skip_serializing_if = "Option::is_none")]
    label: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    tooltip: Option<String>,
}
