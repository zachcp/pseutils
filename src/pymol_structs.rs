//! This module provides data structures and utilities for working with PyMOL session files (PSE).
//!
//! It includes definitions for various PyMOL objects, such as molecules, selectors, and other
//! session-related data. The main structure, `PSEData`, represents the overall content of a
//! PSE file, while other structures like `PyObjectMolecule`, `CoordSet`, `AtomInfo`, and `Bond`
//! represent specific components within the session.
//!
//! This module is designed to facilitate the deserialization of PSE files using the serde and
//! serde_pickle libraries, allowing for easy parsing and manipulation of PyMOL session data
//! in Rust applications.
//!
//!
//! /// PyObject Serialization
//!  https://github.com/schrodinger/pymol-open-source/blob/03d7a7fcf0bd95cd93d710a1268dbace2ed77765/layer1/PyMOLObject.cpp#L681
//!
//! // - Object:
//!      - Gadget
//!      - Molecule
//!      - Dist
//!      - Map
//!      - Mesh
//!      - Slice
//!      - Surface
//!      - CGO
//!      - Alignment
//!      - Group
//!      - Volume
//!      - Callback
//!      - Curve
//!  - Selection
//!

use serde::{Deserialize, Deserializer, Serialize};
use serde_pickle::de::{from_reader, DeOptions};
use serde_pickle::{de, from_value, Value};
use std::io::{self, Read};
use std::{collections::HashMap, fs::File};

#[derive(Debug, Deserialize, Serialize)]
pub struct PSEData {
    pub version: i32,
    main: Vec<i64>,
    colors: Vec<i32>,
    color_ext: Vec<i32>,
    unique_settings: Vec<i32>,
    selector_secrets: Vec<i32>,
    editor: Vec<i32>,
    view: Vec<f32>,
    view_dict: HashMap<String, String>,
    #[serde(with = "serde_bytes")]
    wizard: Vec<u8>,
    moviescenes: Vec<Vec<i32>>,
    settings: Vec<(i32, i32, CustomValue)>,
    // movie: (i32, i32, Vec<f32>, i32, None, None, None),
    movie: (
        i32,
        i32,
        Vec<f32>,
        i32,
        Option<bool>, // this will probably need to be modified
        Option<bool>,
        Option<bool>,
    ),
    // not needed?
    // session: HashMap<String, Value>,
    cache: Vec<usize>,
    // name is the trickiest bit
    names: Vec<Option<SessionName>>,
}

impl PSEData {
    pub fn load(file_path: &str) -> Result<Self, Box<dyn std::error::Error>> {
        let mut file = File::open(file_path)?;
        let mut buffer = Vec::new();
        file.read_to_end(&mut buffer)?;

        let options = DeOptions::new()
            .replace_unresolved_globals()
            .decode_strings();

        // https://github.com/birkenfeld/serde-pickle/issues/12#issuecomment-771974607
        // let v: serde_pickle::Value = serde_pickle::from_reader(file).unwrap();
        // let d: Message = serde_pickle::from_value(v).unwrap();
        let pse_data_vals: serde_pickle::Value = from_reader(&buffer[..], options).unwrap();
        let pse_data: PSEData = serde_pickle::from_value(pse_data_vals).unwrap();
        Ok(pse_data)
    }

    pub fn to_json(&self, file_path: &str) -> Result<(), Box<dyn std::error::Error>> {
        let json = serde_json::to_string_pretty(self)?;
        std::fs::write(file_path, json)?;
        Ok(())
    }
}

#[derive(Debug, Deserialize, Serialize)]
struct SessionName {
    name: String,
    object: i32,
    visible: i32,
    unused: Option<bool>,
    unused2: i32,
    // SelectorAsPyList
    // https://github.com/schrodinger/pymol-open-source/blob/03d7a7fcf0bd95cd93d710a1268dbace2ed77765/layer3/Selector.cpp#L2926
    // list of lists
    // String: Name of the Object
    // Vec1: Index Object ( from VLA list )
    // Vec2: Tag Object ( from VLA list )
    // selector: Vec<(String, Vec<i32>, Vec<i32>)>, // this is there the selection bits are
    data: PymolSessionObjectData,
    // data: PyObjectMolecule,
    group: String,
}

#[derive(Debug, Serialize)]
#[serde(untagged)]
enum PymolSessionObjectData {
    PyObjectMolecule(PyObjectMolecule),
    SessionSelector(SessionSelector),
    // MolVariant(PyObjectMolecule),
    // SessionVariant(SessionSelector),
}

impl<'de> Deserialize<'de> for PymolSessionObjectData {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let value = Value::deserialize(deserializer)?;

        // Debug print
        println!("Deserialized value: {:?}", value);

        // Try to deserialize as PyObjectMolecule
        match from_value::<PyObjectMolecule>(value.clone()) {
            Ok(molecule) => return Ok(PymolSessionObjectData::PyObjectMolecule(molecule)),
            Err(_) => {} // If it fails, we'll try the next option
        }

        // If that fails, try to deserialize as SessionSelector
        match from_value::<SessionSelector>(value.clone()) {
            Ok(selection) => return Ok(PymolSessionObjectData::SessionSelector(selection)),
            Err(_) => {} // If it fails, we'll return an error
        }

        // If both fail, return a generic error
        // Err(String::from("We are unable to serialize this Value"));
        Err(panic!("Problem opening the file"))
    }
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(untagged)]
enum CustomValue {
    Integer(i64),
    Float(f64),
    String(String),
    Boolean(bool),
    // Add more variants if needed, e.g., None for Python's None
}

#[derive(Debug, Deserialize, Serialize)]
struct PyObject {
    object_type: i32,
    name: String,
    color: i32,
    vis_rep: i32,
    extent_min: [f32; 3],
    extent_max: [f32; 3],
    extent_flag: i32,
    ttt_flag: i32,
    setting: Option<bool>, // this is a hack
    enabled: i32,
    render_context: i32,
    ttt: Vec<f32>,
    n_frame: i32,
    view_elem: Option<bool>, //hack
}

/// Coord Set
/// https://github.com/schrodinger/pymol-open-source/blob/master/layer2/CoordSet.cpp#L363
#[derive(Debug, Deserialize, Serialize)]
struct CoordSet {
    n_index: i32,
    n_at_index: i32,
    coord: Vec<f32>,
    idx_to_atm: Vec<i32>,
    atm_to_idx: Option<Vec<i32>>,
    name: String,
    // PyObject *ObjectStateAsPyList(CObjectState * I)
    // https://github.com/schrodinger/pymol-open-source/blob/master/layer1/PyMOLObject.cpp#L1391
    // object_state: Vec<Option<Vec<f64>>>, //
    // *SettingAsPyList
    // https://github.com/schrodinger/pymol-open-source/blob/master/layer1/Setting.cpp#L962
    setting: Vec<Option<bool>>, // punting on this
    lab_pos: Option<bool>,      // might be wrong
    field_9: Option<bool>,      // probably wrong...
    // /// https://github.com/schrodinger/pymol-open-source/blob/master/layer1/CGO.cpp#L220
    sculpt_cgo: Option<(i32, Vec<f32>)>,           //
    atom_state_settings: Option<Vec<Option<i32>>>, //
    // /// https://github.com/schrodinger/pymol-open-source/blob/master/layer1/Symmetry.cpp#L30
    symmetry: Option<Vec<(((i32, i32, i32), (i32, i32, i32)), String)>>, //
}

#[derive(Debug, Deserialize, Serialize)]
struct AtomInfo {
    resv: i32,
    chain: String,
    alt: String,
    resi: String,
    segi: String,
    resn: String,
    name: String,
    elem: String,
    text_type: String,
    label: String,
    ss_type: String,
    is_hydrogen: i8, // this is a boolean
    custom_type: i32,
    priority: i32,
    b: f64,
    q: f64,
    vdw: f64,
    partial_charge: f64,
    formal_charge: i32,
    hetatm: i8, // this is a boolean
    vis_rep: i32,
    color: i32,
    id: i32,
    cartoon: i32,
    flags: i64,
    is_bonded: i8,  // this is a boolean
    chem_flag: i32, // this is a boolean
    geom: i32,
    valence: i32,
    is_masked: i8,    // this is a boolean
    is_protected: i8, // this is a boolean
    protons: i32,
    unique_id: i64,
    stereo: i8,
    discrete_state: i32,
    elec_radius: f64,
    rank: i32,
    hb_donor: i8,    // this is a boolean
    hb_acceptor: i8, // this is a boolean
    atomic_color: i32,
    has_setting: i8, // this is a boolean
    anisou_1: f32,
    anisou_2: f32,
    anisou_3: f32,
    anisou_4: f32,
    anisou_5: f32,
    anisou_6: f32,
    custom: String,
}

#[derive(Debug, Deserialize, Serialize)]
struct Bond {
    index_1: i32,
    index_2: i32,
    order: i32,
    id: i32,
    stereo: i32,
    unique_id: i32,
    has_setting: i32,
    // todo hhandle arrity 7 or arrity 8 with specific symmetry info
    // Symmetry operation of the second atom.
    // symop_2: Option<String>,
}

/// PyObjectMolecule
/// https://github.com/schrodinger/pymol-open-source/blob/master/layer2/ObjectMolecule2.cpp#L3524
/// ObjectMolecule
/// https://github.com/schrodinger/pymol-open-source/blob/03d7a7fcf0bd95cd93d710a1268dbace2ed77765/layer2/ObjectMolecule.h#L58
///
#[derive(Debug, Deserialize, Serialize)]
struct PyObjectMolecule {
    object: PyObject,
    n_cset: i32,
    n_bond: i32,
    n_atom: i32,
    coord_set: Vec<CoordSet>,
    cs_tmpl: Option<Vec<CoordSet>>,
    // /// https://github.com/schrodinger/pymol-open-source/blob/master/layer2/ObjectMolecule2.cpp#L3037
    bond: Vec<Bond>,
    // /// https://github.com/schrodinger/pymol-open-source/blob/master/layer2/ObjectMolecule2.cpp#L3248
    // /// https://github.com/schrodinger/pymol-open-source/blob/master/layer2/AtomInfo.cpp#L792
    atom: Vec<AtomInfo>,
    discrete_flag: i32,
    n_discrete: i32,
    symmetry: Option<(([f32; 3], [f32; 3]), String)>, // crystal space group and name
    cur_cset: i32,
    bond_counter: i32,
    atom_counter: i32,
    discrete_atm_to_idx: Option<Vec<i32>>,
    dcs: Option<Vec<i32>>,
}

#[derive(Debug, Deserialize, Serialize)]
struct SessionSelector {
    name: String,
    data1: Vec<i32>,
    data2: Vec<i32>,
}

// Todo:
//
// PyObjectGadget
// https://github.com/schrodinger/pymol-open-source/blob/master/layer3/Executive.cpp#L5237
// https://github.com/schrodinger/pymol-open-source/blob/master/layer2/ObjectGadget.cpp#L386
// struct PyObjectGadget {}
// struct PyObjectDist {}
// struct PyObjectMap {}
// struct PyObjectMesh {}
// struct PyObjectSurface {}
// struct PyObjectCGO {}
// struct PyObjectAlignment {}
// struct PyObjectGroup {}
// struct PyObjectVolume {}
// struct PyObjectCallback {}
// struct PyObjectCurve {}
