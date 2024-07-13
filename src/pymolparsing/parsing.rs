use pdbtbx::{self, PDB};
use serde::{Deserialize, Deserializer, Serialize};
use serde_pickle::{from_value, Value};

#[derive(Debug, Deserialize, Serialize)]
pub struct SessionName {
    pub name: String,
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
    pub data: PymolSessionObjectData,
    group: String,
}

#[derive(Debug, Serialize)]
#[serde(untagged)]
pub enum PymolSessionObjectData {
    PyObjectMolecule(PyObjectMolecule),
    SessionSelectorList(SessionSelectorList),
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
        // println!("Deserialized value: {:?}", value);

        // Try to deserialize as PyObjectMolecule
        match from_value::<PyObjectMolecule>(value.clone()) {
            Ok(molecule) => return Ok(PymolSessionObjectData::PyObjectMolecule(molecule)),
            Err(_) => {} // If it fails, we'll try the next option
        }

        // println!(
        //     "Did not serialize as a molecule. Not trying as a session: {:?}",
        //     value
        // );

        // If that fails, try to deserialize as SessionSelector
        match from_value::<SessionSelectorList>(value.clone()) {
            Ok(selection) => return Ok(PymolSessionObjectData::SessionSelectorList(selection)),
            Err(_) => {} // If it fails, we'll return an error
        }

        // If both fail, return a generic error
        // Err(String::from("We are unable to serialize this Value"));
        Err(panic!("Problem opening the file"))
    }
}

/// PyObjectMolecule
/// https://github.com/schrodinger/pymol-open-source/blob/master/layer2/ObjectMolecule2.cpp#L3524
/// ObjectMolecule
/// https://github.com/schrodinger/pymol-open-source/blob/03d7a7fcf0bd95cd93d710a1268dbace2ed77765/layer2/ObjectMolecule.h#L58
///
#[derive(Debug, Deserialize, Serialize)]
pub struct PyObjectMolecule {
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
impl PyObjectMolecule {
    pub fn to_pdb(&self) -> PDB {
        // add atoms
        // add bonds
        // add chains
        // add models
        //
        let pdb = PDB::new();

        println!("First Atom: {:?}", &self.atom[1]);
        // let model = Model.new();
        // let chain = Chain.new();
        // let residue = Residue.new();
        // let atom = Atom.new();
        // let bond = Bond.new();
        //
        //
        //
        pdb
    }
}

#[derive(Debug, Deserialize, Serialize)]
struct SessionSelectorList(Vec<SessionSelector>);

#[derive(Debug, Serialize)]
struct SessionSelector {
    id: String,
    values1: Vec<i64>,
    values2: Vec<i64>,
}

// You might need a custom Deserialize implementation for SessionSelector
impl<'de> Deserialize<'de> for SessionSelector {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let (id, values1, values2) = Deserialize::deserialize(deserializer)?;
        Ok(SessionSelector {
            id,
            values1,
            values2,
        })
    }
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(untagged)]
pub enum CustomValue {
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
    pub name: String,
    pub elem: String,
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
    pub hetatm: i8, // this is a boolean
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
impl AtomInfo {
    pub fn is_hetero(&self) -> bool {
        match self.hetatm {
            1 => true,
            0 => false,
            _ => false,
        }
    }
    pub fn to_pdbtbx_atom(&self) -> pdbtbx::Atom {
        let formal_charge = self.formal_charge as isize;

        let atom = pdbtbx::Atom::new(
            self.is_hetero(), // hetero
            0,                // serial_number
            &self.name,       // atom_name
            0.0,              // x Todo
            0.0,              // y Todo
            0.0,              // z Todo
            0.0,              // occupancy? Todo
            self.b,           // b-factor
            &self.elem,       // element
            formal_charge,    // charge: todo: is this the right charge?
        );
        atom.unwrap()
    }
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
