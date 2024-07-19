//! This module provides structures and implementations for handling PyMOL session data.
//!
//! It includes definitions for various PyMOL objects such as molecules, selectors,
//! and atom information. The module also provides functionality for converting
//! PyMOL data structures to PDB format using the `pdbtbx` crate.
//!
//! Key structures:
//! - `SessionName`: Represents a named session object
//! - `PyObjectMolecule`: Represents a molecule in a PyMOL session
//! - `SessionSelector`: Represents a selection in a PyMOL session
//! - `AtomInfo`: Contains detailed information about individual atoms
//!
//! This module is designed to work with serialized PyMOL session data and
//! provides methods for deserializing and manipulating this data.
//!
//! ## Links
//!
//! - [Molecule Experter](https://github.com/schrodinger/pymol-open-source/blob/master/layer3/MoleculeExporter.cpp#L1627)
//! - [PymolMoleculeExporter](https://github.com/schrodinger/pymol-open-source/blob/03d7a7fcf0bd95cd93d710a1268dbace2ed77765/layer4/Cmd.cpp#L3877)
//! - [PDB Exporter](https://github.com/schrodinger/pymol-open-source/blob/master/layer3/MoleculeExporter.cpp#L1627)
//! - [MoleculeExporterPDB](https://github.com/schrodinger/pymol-open-source/blob/master/layer3/MoleculeExporter.cpp#L439)
//! - [Gadget_01](https://github.com/schrodinger/pymol-open-source/blob/master/layer3/Executive.cpp#L5237)
//! - [Gadget_02](https://github.com/schrodinger/pymol-open-source/blob/master/layer2/ObjectGadget.cpp#L386)
//!
//! exporter->init(G);
//! exporter->setMulti(multi);
//! exporter->setRefObject(ref_object, ref_state);
//! exporter->execute(sele, state);
//!
//! CoordSetAtomToPDBStrVLA
//! Variables:
//! - m_iter: SeleCoordIterator m_iter;
//! Atoms:
//! how do we check for multiple objects?
//! m_iter.obj defines the object. By number? by name?
//! - iterate through the coordinates
//! - check for multi
//! update transformation matrices
//! updateMatrix(m_mat_full, true);
//! updateMatrix(m_mat_move, false);
//! beginCoordSet();
//! m_last_cs = m_iter.cs;
//! for bonds
//! if (!m_tmpids[m_iter.getAtm()]) {
//! m_id = m_retain_ids ? m_iter.getAtomInfo()->id : (m_id + 1);
//!  m_tmpids[m_iter.getAtm()] = m_id;
use crate::molviewspec::nodes::{ComponentExpression, ComponentSelector};
use itertools::Itertools;
use pdbtbx::{self, Residue, PDB};
use serde::{Deserialize, Deserializer, Serialize};
use serde_pickle::{from_value, Value};
use serde_repr::{Deserialize_repr, Serialize_repr};

/// AtomInfo
///
/// This struct contains various properties of an atom, including its position,
/// chemical properties, and visualization settings.
///
/// # Fields
///
/// * `resv` - Residue sequence number
/// * `chain` - Chain identifier
/// * `alt` - Alternate location indicator
/// * `resi` - Residue identifier
/// * `segi` - Segment identifier
/// * `resn` - Residue name
/// * `name` - Atom name
/// * `elem` - Element symbol
/// * `text_type` - Text type
/// * `label` - Label text
/// * `ss_type` - Secondary structure type
/// * `is_hydrogen` - Flag indicating if the atom is hydrogen
/// * `custom_type` - Custom type identifier
/// * `priority` - Priority value
/// * `b` - B-factor (temperature factor)
/// * `q` - Occupancy
/// * `vdw` - Van der Waals radius
/// * `partial_charge` - Partial charge
/// * `formal_charge` - Formal charge
/// * `hetatm` - Flag indicating if the atom is a heteroatom
/// * `vis_rep` - Visualization representation
/// * `color` - Color index
/// * `id` - Atom ID
/// * `cartoon` - Cartoon representation type
/// * `flags` - Various flags
/// * `is_bonded` - Flag indicating if the atom is bonded
/// * `chem_flag` - Chemical flag
/// * `geom` - Geometry type
/// * `valence` - Valence
/// * `is_masked` - Flag indicating if the atom is masked
/// * `is_protected` - Flag indicating if the atom is protected
/// * `protons` - Number of protons
/// * `unique_id` - Unique identifier
/// * `stereo` - Stereochemistry indicator
/// * `discrete_state` - Discrete state
/// * `elec_radius` - Electronic radius
/// * `rank` - Rank
/// * `hb_donor` - Hydrogen bond donor flag
/// * `hb_acceptor` - Hydrogen bond acceptor flag
/// * `atomic_color` - Atomic color
/// * `has_setting` - Flag indicating if the atom has custom settings
/// * `anisou_1` to `anisou_6` - Anisotropic temperature factors
/// * `custom` - Custom data string
///
#[derive(Debug, Deserialize, Serialize)]
pub struct AtomInfo {
    pub resv: i32,
    pub chain: String,
    pub alt: String,
    pub resi: String,
    pub segi: String,
    pub resn: String,
    pub name: String,
    pub elem: String,
    pub text_type: String,
    pub label: String,
    pub ss_type: String,
    pub is_hydrogen: i8,
    pub custom_type: i32,
    pub priority: i32,
    pub b: f64,
    pub q: f64,
    pub vdw: f64,
    pub partial_charge: f64,
    pub formal_charge: i32,
    pub hetatm: i8,
    pub vis_rep: i32,
    pub color: i32,
    pub id: i32,
    pub cartoon: i32,
    pub flags: i64,
    pub is_bonded: i8,
    pub chem_flag: i32,
    pub geom: i32,
    pub valence: i32,
    pub is_masked: i8,
    pub is_protected: i8,
    pub protons: i32,
    pub unique_id: i64,
    pub stereo: i8,
    pub discrete_state: i32,
    pub elec_radius: f64,
    pub rank: i32,
    pub hb_donor: i8,
    pub hb_acceptor: i8,
    pub atomic_color: i32,
    pub has_setting: i8,
    pub anisou_1: f32,
    pub anisou_2: f32,
    pub anisou_3: f32,
    pub anisou_4: f32,
    pub anisou_5: f32,
    pub anisou_6: f32,
    pub custom: String,
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

/// Bond Structure
///
/// Represents a chemical bond between two atoms in a molecule.
///
/// # Fields
///
/// * `index_1` - Index of the first atom in the bond
/// * `index_2` - Index of the second atom in the bond
/// * `order` - Bond order (e.g., single, double, triple)
/// * `id` - Unique identifier for the bond
/// * `stereo` - Stereochemistry information for the bond
/// * `unique_id` - Another unique identifier for the bond
/// * `has_setting` - Flag indicating if the bond has custom settings
#[derive(Debug, Deserialize, Serialize)]
pub struct Bond {
    pub index_1: i32,
    pub index_2: i32,
    pub order: i32,
    pub id: i32,
    pub stereo: i32,
    pub unique_id: i32,
    pub has_setting: i32,
    // todo handle arity 7 or arity 8 with specific symmetry info
    // Symmetry operation of the second atom.
    // symop_2: Option<String>,
}

/// Coord Set: Class for storage of coordinates
///
/// Links:
/// - [pymol_coordset](https://github.com/schrodinger/pymol-open-source/blob/master/layer2/CoordSet.cpp#L363)
/// - [pymol_coordset_settings]( https://github.com/schrodinger/pymol-open-source/blob/master/layer1/Setting.cpp#L962)
/// - [CGO](https://github.com/schrodinger/pymol-open-source/blob/master/layer1/CGO.cpp#L220)
/// - [symettry_settings](https://github.com/schrodinger/pymol-open-source/blob/master/layer1/Symmetry.cpp#L30)
///
#[derive(Debug, Deserialize, Serialize)]
pub struct CoordSet {
    pub n_index: i32,         // 1519
    n_at_index: i32,          // 1519
    pub coord: Vec<f32>,      // len -== 4556 ( 1519 *3 )
    pub idx_to_atm: Vec<i32>, // 1 - 1518
    pub atm_to_idx: Option<Vec<i32>>,
    pub name: String,
    pub setting: Vec<Option<bool>>,          // punting on this
    pub lab_pos: Option<bool>,               // might be wrong
    field_9: Option<bool>,                   // probably wrong...
    pub sculpt_cgo: Option<(i32, Vec<f32>)>, //
    pub atom_state_settings: Option<Vec<Option<i32>>>, //
    pub symmetry: Option<Vec<(((i32, i32, i32), (i32, i32, i32)), String)>>,
}

/// Custom Value
///
/// needed for the settings triplet.
///
#[derive(Debug, Serialize, Deserialize)]
#[serde(untagged)]
pub enum CustomValue {
    Integer(i64),
    Float(f64),
    String(String),
    Boolean(bool),
}

/// PyObject
///
/// General settings object
///
#[derive(Debug, Deserialize, Serialize)]
pub struct PyObject {
    pub object_type: i32,
    pub name: String,
    pub color: i32,
    pub vis_rep: i32,
    pub extent_min: [f32; 3],
    pub extent_max: [f32; 3],
    pub extent_flag: i32,
    pub ttt_flag: i32,
    pub setting: Option<bool>, // this is a hack
    pub enabled: i32,
    pub render_context: i32,
    pub ttt: Vec<f32>,
    pub n_frame: i32,
    pub view_elem: Option<bool>, //hack
}

/// PyObjectMolecule: Represents a molecule in PyMOL.
///
/// ## Link
///
/// - [pymol code](https://github.com/schrodinger/pymol-open-source/blob/master/layer2/ObjectMolecule2.cpp#L3524)
/// - [ObjectMolecule](https://github.com/schrodinger/pymol-open-source/blob/03d7a7fcf0bd95cd93d710a1268dbace2ed77765/layer2/ObjectMolecule.h#L58)
/// - [Bond](https://github.com/schrodinger/pymol-open-source/blob/master/layer2/ObjectMolecule2.cpp#L3037)
/// - [AtomInfo](https://github.com/schrodinger/pymol-open-source/blob/master/layer2/ObjectMolecule2.cpp#L3248)
/// - [AtomInfo2](https://github.com/schrodinger/pymol-open-source/blob/master/layer2/AtomInfo.cpp#L792)
///
/// ## Fields
///
/// * `object` - The base PyObject information
/// * `n_cset` - Number of coordinate sets
/// * `n_bond` - Number of bonds
/// * `n_atom` - Number of atoms
/// * `coord_set` - Vector of coordinate sets
/// * `cs_tmpl` - Optional template coordinate set
/// * `bond` - Vector of bonds
/// * `atom` - Vector of atom information
/// * `discrete_flag` - Flag for discrete representation
/// * `n_discrete` - Number of discrete objects
/// * `symmetry` - Optional symmetry information
/// * `cur_cset` - Current coordinate set index
/// * `bond_counter` - Counter for bonds
/// * `atom_counter` - Counter for atoms
/// * `discrete_atm_to_idx` - Optional mapping of discrete atoms to indices
/// * `dcs` - Optional discrete coordinate set information
///
#[derive(Debug, Deserialize, Serialize)]
pub struct PyObjectMolecule {
    pub object: PyObject,
    pub n_cset: i32,
    pub n_bond: i32,
    pub n_atom: i32,
    /// Vector of Coordinates
    pub coord_set: Vec<CoordSet>,
    pub cs_tmpl: Option<Vec<CoordSet>>,
    pub bond: Vec<Bond>,
    pub atom: Vec<AtomInfo>,
    pub discrete_flag: i32,
    pub n_discrete: i32,
    pub symmetry: Option<(([f32; 3], [f32; 3]), String)>, // crystal space group and name
    pub cur_cset: i32,
    pub bond_counter: i32,
    pub atom_counter: i32,
    pub discrete_atm_to_idx: Option<Vec<i32>>,
    pub dcs: Option<Vec<i32>>,
}
impl PyObjectMolecule {
    pub fn get_name(&self) -> String {
        self.object.name.to_string()
    }
    /// Create a PDBTBX::Atom from the pymol object datastructure
    pub fn get_atom(&self, atm_idx: i32) -> pdbtbx::Atom {
        // find atom ids and coordinates in the CoordSet
        // find the remaining atom info in the AtomInfo Vector
        let cset = &self.coord_set;
        // println!("{:?}", cset);
        let atom_coords = &cset[0].coord; // note there may be more than one coord set.... Todo.

        // println!("{:?}", atom_coords);
        // coords are stored in a 1D vector of x,y,z,x,y,x,z,x,y,z
        let base_coord = (3 * atm_idx) as usize;
        // println!("{}", atm_idx);
        // println!("{}", base_coord);
        let x_coord = atom_coords[base_coord];
        // println!("{}", x_coord);
        let y_coord = atom_coords[base_coord + 1];
        // println!("{}", y_coord);
        let z_coord = atom_coords[base_coord + 2];
        // println!("{}", z_coord);
        // println!("{}, {}, {}", x_coord, y_coord, z_coord);

        let atom_info = &self.atom.iter().find(|atm| atm.id == atm_idx + 1).unwrap(); // note that the atom in the atom vector seem to be 1-indexed.
        let formal_charge = atom_info.formal_charge as isize;
        let serial_number = atom_info.id as usize;

        let atom = pdbtbx::Atom::new(
            atom_info.is_hetero(),  // hetero
            serial_number,          // serial_number: Note: I am not sure this is correct just yet.
            atom_info.name.clone(), // atom_name
            x_coord.into(),         // x
            y_coord.into(),         // y
            z_coord.into(),         // z
            0.0,                    // occupancy? Todo
            atom_info.b,            // b-factor
            atom_info.elem.clone(), // element
            formal_charge,          // charge: todo: is this the right charge?
        );
        atom.unwrap()
    }
    /// Get unique chain names
    pub fn get_chains(&self) -> Vec<String> {
        self.atom
            .iter()
            .filter_map(|atm| Some(atm.chain.clone()))
            .unique() // from itertools
            .collect()
    }
    /// Get each residue by chain.
    pub fn get_residues_by_chain(&self, chain: String) -> Vec<i32> {
        self.atom
            .iter()
            .filter(|atm| atm.chain == chain)
            .map(|atm| atm.resv)
            .unique()
            .collect()
    }
    pub fn get_unit_cell_symmetry(&self) -> (pdbtbx::UnitCell, pdbtbx::Symmetry) {
        let symmetry = &self.symmetry.clone().expect("Expected a symmetry group.");
        let (([a, b, c], [alpha, beta, gamma]), sym_group) = symmetry;
        println!(
            "{}, {}, {}, {}, {}, {}, {}",
            a, b, c, alpha, beta, gamma, sym_group
        );

        // let unitcell = pdbtbx::UnitCell::new(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        let unitcell = pdbtbx::UnitCell::new(
            *a as f64,
            *b as f64,
            *c as f64,
            *alpha as f64,
            *beta as f64,
            *gamma as f64,
        );

        let pdbsym = pdbtbx::Symmetry::new(sym_group).expect("Invalid Symmetry group");

        (unitcell, pdbsym)
    }
    /// Get each residue by chain.
    pub fn create_residue(&self, chain: String, residue_number: i32) -> pdbtbx::Residue {
        let atoms: Vec<&AtomInfo> = self
            .atom
            .iter()
            .filter(|atm| atm.chain == chain && atm.resv == residue_number)
            .collect();

        let resv = residue_number as isize;
        let res_name = atoms[0].name.clone();
        let mut residue = pdbtbx::Residue::new(resv, None, None).expect("Couldn't create residue");
        let mut conformer =
            pdbtbx::Conformer::new(res_name, None, None).expect("Couldn't create Conformer");

        for atom in atoms {
            // coordinate vector is zero-indexed so we need to subtract 1
            let atom = &self.get_atom(atom.id - 1);
            conformer.add_atom(atom.clone());
        }

        residue.add_conformer(conformer);
        residue
    }
    pub fn create_chain(&self, chain: String) -> pdbtbx::Chain {
        let mut new_chain = pdbtbx::Chain::new(chain.clone()).unwrap();

        let residues: Vec<Residue> = self
            .get_residues_by_chain(chain.clone())
            .iter()
            .map(|res_num| self.create_residue(chain.clone(), *res_num))
            .collect();

        // index out of bounds: the len is 4557 but the index is 4557
        for res in residues {
            new_chain.add_residue(res.clone())
        }

        new_chain
    }
    // Create a pdbtbx::PDB
    pub fn to_pdb(&self) -> PDB {
        // Create a Model. Need to fix this later if theres multiple models
        let mut model = pdbtbx::Model::new(1);
        let chains: Vec<pdbtbx::Chain> = self
            .get_chains()
            .iter()
            .map(|chainid| self.create_chain(chainid.to_string()))
            .collect();

        for chain in chains {
            model.add_chain(chain);
        }

        // Create PDB from Models
        let mut pdb = PDB::new();
        pdb.add_model(model);

        // Add Bonds Here
        for bond in &self.bond {
            // Todo: proper pymol bond--> pdbtbx bond
            pdb.add_bond(
                (bond.index_1 as usize, None),
                (bond.index_2 as usize, None),
                pdbtbx::Bond::Covalent,
            );
        }

        // Add Name/ Identifier
        let identifier = self.get_name().clone();
        pdb.identifier = Some(identifier);

        // Add Unit Cell and Symmetrey Info
        let (unit, sym) = self.get_unit_cell_symmetry();
        pdb.unit_cell = Some(unit);
        pdb.symmetry = Some(sym);

        pdb
    }
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

#[derive(Debug, Deserialize, Serialize)]
pub struct SessionName {
    pub name: String,
    pub object: i32,
    pub visible: i32,
    unused: Option<bool>,
    unused2: i32,
    pub data: PymolSessionObjectData,
    pub group: String,
}

#[derive(Debug, Serialize)]
pub struct SessionSelector {
    // SelectorAsPyList
    // https://github.com/schrodinger/pymol-open-source/blob/03d7a7fcf0bd95cd93d710a1268dbace2ed77765/layer3/Selector.cpp#L2926
    // list of lists
    // String: Name of the Object
    // Vec1: Index Object ( from VLA list )
    // Vec2: Tag Object ( from VLA list )
    // selector: Vec<(String, Vec<i32>, Vec<i32>)>, // this is there the selection bits are
    pub id: String,
    pub atom_index: Vec<i64>,
    pub atom_tag: Vec<i64>,
}
impl SessionSelector {
    pub fn to_component(&self) -> ComponentSelector {
        let mut expression_list: Vec<ComponentExpression> = vec![];
        for idx in &self.atom_index {
            let idx32: i32 = *idx as i32;
            expression_list.push(ComponentExpression {
                // internal representation of selection is available as the atom_index
                atom_index: Some(idx32),
                ..Default::default()
            });
        }
        ComponentSelector::ExpressionList(expression_list)
    }
}

// You might need a custom Deserialize implementation for SessionSelector
impl<'de> Deserialize<'de> for SessionSelector {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let (id, atom_index, atom_tag) = Deserialize::deserialize(deserializer)?;
        Ok(SessionSelector {
            id,
            atom_index,
            atom_tag,
        })
    }
}

#[derive(Debug, Deserialize, Serialize)]
pub struct SessionSelectorList(Vec<SessionSelector>);
impl SessionSelectorList {
    pub fn get_selectors(&self) -> &Vec<SessionSelector> {
        &self.0
    }
}

#[derive(Serialize_repr, Deserialize_repr, PartialEq, Debug)]
#[repr(u32)]
pub enum SettingsEnum {
    BondingVdwCutoff = 0,
    MinMeshSpacing = 1,
    DotDensity = 2,
    DotMode = 3,
    SolventRadius = 4,
    SelCounter = 5,
    BgRgb = 6,
    Ambient = 7,
    Direct = 8,
    Reflect = 9,
    Light = 10,
    Power = 11,
    Antialias = 12,
    CavityCull = 13,
    GlAmbient = 14,
    SingleImage = 15,
    MovieDelay = 16,
    RibbonPower = 17,
    RibbonPowerB = 18,
    RibbonSampling = 19,
    RibbonRadius = 20,
    StickRadius = 21,
    HashMax = 22,
    Orthoscopic = 23,
    SpecReflect = 24,
    SpecPower = 25,
    SweepAngle = 26,
    SweepSpeed = 27,
    DotHydrogens = 28,
    DotRadius = 29,
    RayTraceFrames = 30,
    CacheFrames = 31,
    TrimDots = 32,
    CullSpheres = 33,
    Test1 = 34,
    Test2 = 35,
    SurfaceBest = 36,
    SurfaceNormal = 37,
    SurfaceQuality = 38,
    SurfaceProximity = 39,
    NormalWorkaround = 40,
    StereoAngle = 41,
    StereoShift = 42,
    LineSmooth = 43,
    LineWidth = 44,
    HalfBonds = 45,
    StickQuality = 46,
    StickOverlap = 47,
    StickNub = 48,
    AllStates = 49,
    Pickable = 50,
    AutoShowLines = 51,
    IdleDelay = 52,
    NoIdle = 53,
    FastIdle = 54,
    SlowIdle = 55,
    RockDelay = 56,
    DistCounter = 57,
    DashLength = 58,
    DashGap = 59,
    AutoZoom = 60,
    Overlay = 61,
    Text = 62,
    ButtonMode = 63,
    Valence = 64,
    NonbondedSize = 65,
    LabelColor = 66,
    RayTraceFog = 67,
    SpheroidScale = 68,
    RayTraceFogStart = 69,
    SpheroidSmooth = 70,
    SpheroidFill = 71,
    AutoShowNonbonded = 72,
    CacheDisplay = 73,
    MeshRadius = 74,
    BackfaceCull = 75,
    Gamma = 76,
    DotWidth = 77,
    AutoShowSelections = 78,
    AutoHideSelections = 79,
    SelectionWidth = 80,
    SelectionOverlay = 81,
    StaticSingletons = 82,
    // Unused83 = 83,
    DepthCue = 84,
    Specular = 85,
    Shininess = 86,
    SphereQuality = 87,
    Fog = 88,
    IsomeshAutoState = 89,
    MeshWidth = 90,
    CartoonSampling = 91,
    CartoonLoopRadius = 92,
    CartoonLoopQuality = 93,
    CartoonPower = 94,
    CartoonPowerB = 95,
    CartoonRectLength = 96,
    CartoonRectWidth = 97,
    InternalGuiWidth = 98,
    InternalGui = 99,
    CartoonOvalLength = 100,
    CartoonOvalWidth = 101,
    CartoonOvalQuality = 102,
    CartoonTubeRadius = 103,
    CartoonTubeQuality = 104,
    CartoonDebug = 105,
    RibbonWidth = 106,
    DashWidth = 107,
    DashRadius = 108,
    CgoRayWidthScale = 109,
    LineRadius = 110,
    CartoonRoundHelices = 111,
    CartoonRefineNormals = 112,
    CartoonFlatSheets = 113,
    CartoonSmoothLoops = 114,
    CartoonDumbbellLength = 115,
    CartoonDumbbellWidth = 116,
    CartoonDumbbellRadius = 117,
    CartoonFancyHelices = 118,
    CartoonFancySheets = 119,
    IgnorePdbSegi = 120,
    RibbonThrow = 121,
    CartoonThrow = 122,
    CartoonRefine = 123,
    CartoonRefineTips = 124,
    CartoonDiscreteColors = 125,
    NormalizeCcp4Maps = 126,
    SurfacePoor = 127,
    InternalFeedback = 128,
    CgoLineWidth = 129,
    CgoLineRadius = 130,
    Logging = 131,
    RobustLogs = 132,
    LogBoxSelections = 133,
    LogConformations = 134,
    ValenceSize = 135,
    SurfaceMiserable = 136,
    RayOpaqueBackground = 137,
    Transparency = 138,
    RayTexture = 139,
    RayTextureSettings = 140,
    SuspendUpdates = 141,
    FullScreen = 142,
    SurfaceMode = 143,
    SurfaceColor = 144,
    MeshMode = 145,
    MeshColor = 146,
    AutoIndicateFlags = 147,
    SurfaceDebug = 148,
    RayImproveShadows = 149,
    SmoothColorTriangle = 150,
    RayDefaultRenderer = 151,
    FieldOfView = 152,
    ReflectPower = 153,
    PreserveChempyIds = 154,
    SphereScale = 155,
    TwoSidedLighting = 156,
    SecondaryStructure = 157,
    AutoRemoveHydrogens = 158,
    RaiseExceptions = 159,
    StopOnExceptions = 160,
    Sculpting = 161,
    AutoSculpt = 162,
    SculptVdwScale = 163,
    SculptVdwScale14 = 164,
    SculptVdwWeight = 165,
    SculptVdwWeight14 = 166,
    SculptBondWeight = 167,
    SculptAnglWeight = 168,
    SculptPyraWeight = 169,
    SculptPlanWeight = 170,
    SculptingCycles = 171,
    SphereTransparency = 172,
    SphereColor = 173,
    SculptFieldMask = 174,
    SculptHbOverlap = 175,
    SculptHbOverlapBase = 176,
    LegacyVdwRadii = 177,
    SculptMemory = 178,
    ConnectMode = 179,
    CartoonCylindricalHelices = 180,
    CartoonHelixRadius = 181,
    ConnectCutoff = 182,
    SavePdbSs = 183,
    SculptLineWeight = 184,
    FitIterations = 185,
    FitTolerance = 186,
    BatchPrefix = 187,
    StereoMode = 188,
    CgoSphereQuality = 189,
    PdbLiteralNames = 190,
    WrapOutput = 191,
    FogStart = 192,
    State = 193,
    Frame = 194,
    RayShadow = 195,
    RibbonTraceAtoms = 196,
    Security = 197,
    StickTransparency = 198,
    RayTransparencyShadows = 199,
    SessionVersionCheck = 200,
    RayTransparencySpecular = 201,
    StereDoublePumpMono = 202,
    SphereSolvent = 203,
    MeshQuality = 204,
    MeshSolvent = 205,
    DotSolvent = 206,
    RayShadowFudge = 207,
    RayTriangleFudge = 208,
    DebugPick = 209,
    DotColor = 210,
    MouseLimit = 211,
    MouseScale = 212,
    TransparencyMode = 213,
    ClampColors = 214,
    PymolSpaceMaxRed = 215,
    PymolSpaceMaxGreen = 216,
    PymolSpaceMaxBlue = 217,
    PymolSpaceMinFactor = 218,
    RovingOrigin = 219,
    RovingLines = 220,
    RovingSticks = 221,
    RovingSpheres = 222,
    RovingLabels = 223,
    RovingDelay = 224,
    RovingSelection = 225,
    RovingByres = 226,
    RovingRibbon = 227,
    RovingCartoon = 228,
    RovingPolarContacts = 229,
    RovingPolarCutoff = 230,
    RovingNonbonded = 231,
    FloatLabels = 232,
    RovingDetail = 233,
    RovingNbSpheres = 234,
    RibbonColor = 235,
    CartoonColor = 236,
    RibbonSmooth = 237,
    AutoColor = 238,
    AutoColorNext = 239,
    RayInteriorColor = 240,
    CartoonHighlightColor = 241,
    CoulombUnitsFactor = 242,
    CoulombDielectric = 243,
    RayInteriorShadows = 244,
    RayInteriorTexture = 245,
    RovingMap1Name = 246,
    RovingMap2Name = 247,
    RovingMap3Name = 248,
    RovingMap1Level = 249,
    RovingMap2Level = 250,
    RovingMap3Level = 251,
    RovingIsomesh = 252,
    RovingIsosurface = 253,
    ScenesChanged = 254,
    GaussianBAdjust = 255,
    PdbStandardOrder = 256,
    CartoonSmoothFirst = 257,
    CartoonSmoothLast = 258,
    CartoonSmoothCycles = 259,
    CartoonFlatCycles = 260,
    MaxThreads = 261,
    ShowProgress = 262,
    UseDisplayLists = 263,
    CacheMemory = 264,
    SimplifyDisplayLists = 265,
    RetainOrder = 266,
    PdbHetatmSort = 267,
    PdbUseTerRecords = 268,
    CartoonTraceAtoms = 269,
    RayOversampleCutoff = 270,
    GaussianResolution = 271,
    GaussianBFloor = 272,
    SculptNbInterval = 273,
    SculptTorsWeight = 274,
    SculptTorsTolerance = 275,
    StickBall = 276,
    StickBallRatio = 277,
    StickFixedRadius = 278,
    CartoonTransparency = 279,
    DashRoundEnds = 280,
    HBondMaxAngle = 281,
    HBondCutoffCenter = 282,
    HBondCutoffEdge = 283,
    HBondPowerA = 284,
    HBondPowerB = 285,
    HBondCone = 286,
    SsHelixPsiTarget = 287,
    SsHelixPsiInclude = 288,
    SsHelixPsiExclude = 289,
    SsHelixPhiTarget = 290,
    SsHelixPhiInclude = 291,
    SsHelixPhiExclude = 292,
    SsStrandPsiTarget = 293,
    SsStrandPsiInclude = 294,
    SsStrandPsiExclude = 295,
    SsStrandPhiTarget = 296,
    SsStrandPhiInclude = 297,
    SsStrandPhiExclude = 298,
    MovieLoop = 299,
    PdbRetainIds = 300,
    PdbNoEndRecord = 301,
    CgoDotWidth = 302,
    CgoDotRadius = 303,
    DeferUpdates = 304,
    NormalizeOMaps = 305,
    SwapDsn6Bytes = 306,
    PdbInsertionsGoFirst = 307,
    RovingOriginZ = 308,
    RovingOriginZCushion = 309,
    SpecularIntensity = 310,
    OverlayLines = 311,
    RayTransparencySpecCut = 312,
    InternalPrompt = 313,
    NormalizeGrdMaps = 314,
    RayBlendColors = 315,
    RayBlendRed = 316,
    RayBlendGreen = 317,
    RayBlendBlue = 318,
    PngScreenGamma = 319,
    PngFileGamma = 320,
    EditorLabelFragments = 321,
    InternalGuiControlSize = 322,
    AutoDss = 323,
    TransparencyPickingMode = 324,
    VirtualTrackball = 325,
    PdbReformatNamesMode = 326,
    RayPixelScale = 327,
    LabelFontId = 328,
    PdbConectAll = 329,
    ButtonModeName = 330,
    SurfaceType = 331,
    DotNormals = 332,
    SessionMigration = 333,
    MeshNormals = 334,
    MeshType = 335,
    DotLighting = 336,
    MeshLighting = 337,
    SurfaceSolvent = 338,
    TriangleMaxPasses = 339,
    RayInteriorReflect = 340,
    InternalGuiMode = 341,
    SurfaceCarveSelection = 342,
    SurfaceCarveState = 343,
    SurfaceCarveCutoff = 344,
    SurfaceClearSelection = 345,
    SurfaceClearState = 346,
    SurfaceClearCutoff = 347,
    SurfaceTrimCutoff = 348,
    SurfaceTrimFactor = 349,
    RayMaxPasses = 350,
    ActiveSelections = 351,
    RayTransparencyContrast = 352,
    SeqView = 353,
    MouseSelectionMode = 354,
    SeqViewLabelSpacing = 355,
    SeqViewLabelStart = 356,
    SeqViewFormat = 357,
    SeqViewLocation = 358,
    SeqViewOverlay = 359,
    AutoClassifyAtoms = 360,
    CartoonNucleicAcidMode = 361,
    SeqViewColor = 362,
    SeqViewLabelMode = 363,
    SurfaceRampAboveMode = 364,
    Stereo = 365,
    WizardPromptMode = 366,
    CoulombCutoff = 367,
    SliceTrackCamera = 368,
    SliceHeightScale = 369,
    SliceHeightMap = 370,
    SliceGrid = 371,
    SliceDynamicGrid = 372,
    SliceDynamicGridResolution = 373,
    PdbInsureOrthogonal = 374,
    RayDirectShade = 375,
    StickColor = 376,
    CartoonPuttyRadius = 377,
    CartoonPuttyQuality = 378,
    CartoonPuttyScaleMin = 379,
    CartoonPuttyScaleMax = 380,
    CartoonPuttyScalePower = 381,
    CartoonPuttyRange = 382,
    CartoonSideChainHelper = 383,
    SurfaceOptimizeSubsets = 384,
    Multiplex = 385,
    TextureFonts = 386,
    PqrNoChainId = 387,
    Animation = 388,
    AnimationDuration = 389,
    SceneAnimation = 390,
    LineStickHelper = 391,
    RayOrthoscopic = 392,
    RibbonSideChainHelper = 393,
    SelectionWidthMax = 394,
    SelectionWidthScale = 395,
    SceneCurrentName = 396,
    Presentation = 397,
    PresentationMode = 398,
    PdbTruncateResidueName = 399,
    SceneLoop = 400,
    SweepMode = 401,
    SweepPhase = 402,
    SceneRestartMovieDelay = 403,
    MouseRestartMovieDelay = 404,
    AngleSize = 405,
    AngleLabelPosition = 406,
    DihedralSize = 407,
    DihedralLabelPosition = 408,
    DeferBuildsMode = 409,
    SeqViewDiscreteByState = 410,
    SceneAnimationDuration = 411,
    Wildcard = 412,
    AtomNameWildcard = 413,
    IgnoreCase = 414,
    PresentationAutoQuit = 415,
    EditorAutoDihedral = 416,
    PresentationAutoStart = 417,
    ValidateObjectNames = 418,
    UnusedBooleanDefTrue = 419,
    AutoShowSpheres = 420,
    SphereMode = 421,
    SpherePointMaxSize = 422,
    SpherePointSize = 423,
    PdbHonorModelNumber = 424,
    RankAssistedSorts = 425,
    RibbonNucleicAcidMode = 426,
    CartoonRingMode = 427,
    CartoonRingWidth = 428,
    CartoonRingColor = 429,
    CartoonRingFinder = 430,
    CartoonTubeCap = 431,
    CartoonLoopCap = 432,
    NvidiaBugs = 433,
    ImageDotsPerInch = 434,
    OpaqueBackground = 435,
    DrawFrames = 436,
    ShowAlphaChecker = 437,
    MatrixMode = 438,
    EditorAutoOrigin = 439,
    SessionFile = 440,
    CgoTransparency = 441,
    LegacyMouseZoom = 442,
    AutoNumberSelections = 443,
    SculptVdwVisMode = 444,
    SculptVdwVisMin = 445,
    SculptVdwVisMid = 446,
    SculptVdwVisMax = 447,
    CartoonLadderMode = 448,
    CartoonLadderRadius = 449,
    CartoonLadderColor = 450,
    CartoonNucleicAcidColor = 451,
    CartoonRingTransparency = 452,
    LabelSize = 453,
    SpecDirect = 454,
    LightCount = 455,
    Light2 = 456,
    Light3 = 457,
    HideUnderscoreNames = 458,
    SelectionRoundPoints = 459,
    DistanceExclusion = 460,
    HBondExclusion = 461,
    LabelShadowMode = 462,
    Light4 = 463,
    Light5 = 464,
    Light6 = 465,
    Light7 = 466,
    LabelOutlineColor = 467,
    RayTraceMode = 468,
    RayTraceGain = 469,
    SelectionVisibleOnly = 470,
    LabelPosition = 471,
    RayTraceDepthFactor = 472,
    RayTraceSlopeFactor = 473,
    RayTraceDiscoFactor = 474,
    RayShadowDecayFactor = 475,
    RayInteriorMode = 476,
    RayLegacyLighting = 477,
    SculptAutoCenter = 478,
    PdbDiscreteChains = 479,
    PdbUnbondCations = 480,
    SculptTriScale = 481,
    SculptTriWeight = 482,
    SculptTriMin = 483,
    SculptTriMax = 484,
    SculptTriMode = 485,
    PdbEchoTags = 486,
    ConnectBonded = 487,
    SpecDirectPower = 488,
    Light8 = 489,
    Light9 = 490,
    RayShadowDecayRange = 491,
    SpecCount = 492,
    SculptMinScale = 493,
    SculptMinWeight = 494,
    SculptMinMin = 495,
    SculptMinMax = 496,
    SculptMaxScale = 497,
    SculptMaxWeight = 498,
    SculptMaxMin = 499,
    SculptMaxMax = 500,
    SurfaceCircumscribe = 501,
    SculptAvdWeight = 502,
    SculptAvdGap = 503,
    SculptAvdRange = 504,
    SculptAvdExcl = 505,
    AsyncBuilds = 506,
    FetchPath = 507,
    CartoonRingRadius = 508,
    RayColorRamps = 509,
    RayHintCamera = 510,
    RayHintShadow = 511,
    StickValenceScale = 512,
    SeqViewAlignment = 513,
    SeqViewUnalignedMode = 514,
    SeqViewUnalignedColor = 515,
    SeqViewFillChar = 516,
    SeqViewFillColor = 517,
    SeqViewLabelColor = 518,
    SurfaceCarveNormalCutoff = 519,
    TraceAtomsMode = 520,
    SessionChanged = 521,
    RayClipShadows = 522,
    MouseWheelScale = 523,
    NonbondedTransparency = 524,
    RaySpecLocal = 525,
    LineColor = 526,
    RayLabelSpecular = 527,
    MeshSkip = 528,
    LabelDigits = 529,
    LabelDistanceDigits = 530,
    LabelAngleDigits = 531,
    LabelDihedralDigits = 532,
    SurfaceNegativeVisible = 533,
    SurfaceNegativeColor = 534,
    MeshNegativeVisible = 535,
    MeshNegativeColor = 536,
    GroupAutoMode = 537,
    GroupFullMemberNames = 538,
    GradientMaxLength = 539,
    GradientMinLength = 540,
    GradientMinSlope = 541,
    GradientNormalMinDot = 542,
    GradientStepSize = 543,
    GradientSpacing = 544,
    GradientSymmetry = 545,
    RayTraceColor = 546,
    GroupArrowPrefix = 547,
    SuppressHidden = 548,
    SessionCompression = 549,
    MovieFps = 550,
    RayTransparencyOblique = 551,
    RayTraceTransCutoff = 552,
    RayTracePersistCutoff = 553,
    RayTransparencyObliquePower = 554,
    RayScatter = 555,
    HBondFromProton = 556,
    AutoCopyImages = 557,
    MoeSeparateChains = 558,
    TransparencyGlobalSort = 559,
    HideLongBonds = 560,
    AutoRenameDuplicateObjects = 561,
    PdbHetatmGuessValences = 562,
    EllipsoidQuality = 563,
    CgoEllipsoidQuality = 564,
    MovieAnimateByFrame = 565,
    RampBlendNearbyColors = 566,
    AutoDeferBuilds = 567,
    EllipsoidProbability = 568,
    EllipsoidScale = 569,
    EllipsoidColor = 570,
    EllipsoidTransparency = 571,
    MovieRock = 572,
    CacheMode = 573,
    DashColor = 574,
    AngleColor = 575,
    DihedralColor = 576,
    GridMode = 577,
    CacheMax = 578,
    GridSlot = 579,
    GridMax = 580,
    CartoonPuttyTransform = 581,
    Rock = 582,
    ConeQuality = 583,
    PdbFormalCharges = 584,
    AtiBugs = 585,
    GeometryExportMode = 586,
    MouseGrid = 587,
    MeshCutoff = 588,
    MeshCarveSelection = 589,
    MeshCarveState = 590,
    MeshCarveCutoff = 591,
    MeshClearSelection = 592,
    MeshClearState = 593,
    MeshClearCutoff = 594,
    MeshGridMax = 595,
    SessionCacheOptimize = 596,
    SdofDragScale = 597,
    SceneButtonsMode = 598,
    SceneButtons = 599,
    MapAutoExpandSym = 600,
    ImageCopyAlways = 601,
    MaxUps = 602,
    AutoOverlay = 603,
    StickBallColor = 604,
    StickHScale = 605,
    SculptPyraInvWeight = 606,
    KeepAlive = 607,
    FitKabsch = 608,
    StereoDynamicStrength = 609,
    DynamicWidth = 610,
    DynamicWidthFactor = 611,
    DynamicWidthMin = 612,
    DynamicWidthMax = 613,
    DrawMode = 614,
    CleanElectroMode = 615,
    ValenceMode = 616,
    ShowFrameRate = 617,
    MoviePanel = 618,
    MouseZScale = 619,
    MovieAutoStore = 620,
    MovieAutoInterpolate = 621,
    MoviePanelRowHeight = 622,
    SceneFrameMode = 623,
    SurfaceCavityMode = 624,
    SurfaceCavityRadius = 625,
    SurfaceCavityCutoff = 626,
    MotionPower = 627,
    MotionBias = 628,
    MotionSimple = 629,
    MotionLinear = 630,
    MotionHand = 631,
    PdbIgnoreConect = 632,
    EditorBondCycleMode = 633,
    MovieQuality = 634,
    LabelAnchor = 635,
    FetchHost = 636,
    DynamicMeasures = 637,
    NeighborCutoff = 638,
    HeavyNeighborCutoff = 639,
    PolarNeighborCutoff = 640,
    SurfaceResidueCutoff = 641,
    SurfaceUseShader = 642,
    CartoonUseShader = 643,
    StickUseShader = 644,
    LineUseShader = 645,
    SphereUseShader = 646,
    UseShaders = 647,
    ShadersFromDisk = 648,
    VolumeBitDepth = 649,
    VolumeColor = 650,
    VolumeLayers = 651,
    VolumeDataRange = 652,
    AutoDeferAtomCount = 653,
    DefaultRefmacNames = 654,
    DefaultPhenixNames = 655,
    DefaultPhenixNoFillNames = 656,
    DefaultBusterNames = 657,
    DefaultFofcMapRep = 658,
    Default2fofcMapRep = 659,
    AtomTypeFormat = 660,
    AutocloseDialogs = 661,
    BgGradient = 662,
    BgRgbTop = 663,
    BgRgbBottom = 664,
    RayVolume = 665,
    RibbonTransparency = 666,
    StateCounterMode = 667,
    CgoUseShader = 668,
    CgoShaderUbColor = 669,
    CgoShaderUbNormal = 670,
    CgoLighting = 671,
    MeshUseShader = 672,
    StickDebug = 673,
    CgoDebug = 674,
    StickRoundNub = 675,
    StickGoodGeometry = 676,
    StickAsCylinders = 677,
    MeshAsCylinders = 678,
    LineAsCylinders = 679,
    RibbonAsCylinders = 680,
    RibbonUseShader = 681,
    ExclDisplayListsShaders = 682,
    DashUseShader = 683,
    DashAsCylinders = 684,
    NonbondedUseShader = 685,
    NonbondedAsCylinders = 686,
    CylindersShaderFilterFaces = 687,
    NbSpheresSize = 688,
    NbSpheresQuality = 689,
    NbSpheresUseShader = 690,
    RenderAsCylinders = 691,
    AlignmentAsCylinders = 692,
    CartoonNucleicAcidAsCylinders = 693,
    CgoShaderUbFlags = 694,
    AntialiasShader = 695,
    OffscreenRenderingMultiplier = 696,
    CylinderShaderFfWorkaround = 697,
    SurfaceColorSmoothing = 698,
    SurfaceColorSmoothingThreshold = 699,
    DotUseShader = 700,
    DotAsSheres = 701,
    AmbientOcclusionMode = 702,
    AmbientOcclusionScale = 703,
    AmbientOcclusionSmooth = 704,
    SmoothHalfBonds = 705,
    AnaglyphMode = 706,
    EditLight = 707,
    SuspendUndo = 708,
    SuspendUndoAtomCount = 709,
    SuspendDeferred = 710,
    PickSurface = 711,
    BgImageFilename = 712,
    BgImageMode = 713,
    BgImageTilesize = 714,
    BgImageLinear = 715,
    LoadObjectPropsDefault = 716,
    LoadAtomPropsDefault = 717,
    LabelPlacementOffset = 718,
    PdbConectNodup = 719,
    LabelConnector = 720,
    LabelConnectorMode = 721,
    LabelConnectorColor = 722,
    LabelConnectorWidth = 723,
    LabelConnectorExtLength = 724,
    LabelBgColor = 725,
    UseGeometryShaders = 726,
    LabelRelativeMode = 727,
    LabelScreenPoint = 728,
    LabelMultilineSpacing = 729,
    LabelMultilineJustification = 730,
    LabelPadding = 731,
    LabelBgTransparency = 732,
    LabelBgOutline = 733,
    RayLabelConnectorFlat = 734,
    DashTransparency = 735,
    PickLabels = 736,
    LabelZTarget = 737,
    SessionEmbedsData = 738,
    VolumeMode = 739,
    Trilines = 740,
    ColladaExportLighting = 741,
    ColladaGeometryMode = 742,
    PrecomputedLighting = 743,
    Chromadepth = 744,
    PseExportVersion = 745,
    CifUseAuth = 746,
    Assembly = 747,
    CifKeepinmemory = 748,
    PseBinaryDump = 749,
    CartoonGapCutoff = 750,
    IgnoreCaseChain = 751,
    ValenceZeroScale = 752,
    ValenceZeroMode = 753,
    AutoShowClassified = 754,
    ColladaBackgroundBox = 755,
    Pick32bit = 756,
    CartoonAllAlt = 757,
    DisplayScaleFactor = 758,
    PickShading = 759,
    FetchTypeDefault = 760,
    EditorAutoMeasure = 761,
    SurfaceSmoothEdges = 762,
    ChemCompCartnUse = 763,
    ColoredFeedback = 764,
    SdfWriteZeroOrderBonds = 765,
    CifMetalcAsZeroOrderBonds = 766,
    SeqViewGapMode = 767,
    InternalGuiNameColorMode = 768,
    OpenvrGuiFov = 769,
    OpenvrGuiAlpha = 770,
    OpenvrGuiUseAlpha = 771,
    OpenvrGuiSceneColor = 772,
    OpenvrGuiSceneAlpha = 773,
    OpenvrGuiBackColor = 774,
    OpenvrGuiBackAlpha = 775,
    OpenvrGuiUseBackdrop = 776,
    OpenvrGuiOverlay = 777,
    OpenvrGuiText = 778,
    OpenvrDisableClipping = 779,
    OpenvrNearPlane = 780,
    OpenvrFarPlane = 781,
    OpenvrCutLaser = 782,
    OpenvrLaserWidth = 783,
    OpenvrGuiDistance = 784,
    CartoonSmoothCylinderCycles = 785,
    CartoonSmoothCylinderWindow = 786,
    IsosurfaceAlgorithm = 787,
    CellCentered = 788,
    HalogenBondDistance = 789,
    HalogenBondAsDonorMinDonorAngle = 790,
    HalogenBondAsDonorMinAcceptorAngle = 791,
    HalogenBondAsAcceptorMinDonorAngle = 792,
    HalogenBondAsAcceptorMinAcceptorAngle = 793,
    HalogenBondAsAcceptorMaxAcceptorAngle = 794,
    SaltBridgeDistance = 795,
    UseTessellationShaders = 796,
    CellColor = 797,
}
// Implement From<i32> for Settings
// impl From<i32> for Settings {
//     fn from(value: i32) -> Self {
//         match value {
//             0 => Settings::BondingVdwCutoff,
//             1 => Settings::MinMeshSpacing,
//             2 => Settings::DotDensity,
//             3 => Settings::DotMode,
//             4 => Settings::SolventRadius,
//             5 => Settings::SelCounter,
//             6 => Settings::BgRgb,
//             7 => Settings::Ambient,
//             8 => Settings::Direct,
//             9 => Settings::Reflect,
//             10 => Settings::Light,
//             11 => Settings::Power,
//             12 => Settings::Antialias,
//             13 => Settings::CavityCull,
//             14 => Settings::GlAmbient,
//             15 => Settings::SingleImage,
//             16 => Settings::MovieDelay,
//             17 => Settings::RibbonPower,
//             18 => Settings::RibbonPowerB,
//             19 => Settings::RibbonSampling,
//             20 => Settings::RibbonRadius,
//             21 => Settings::StickRadius,
//             22 => Settings::HashMax,
//             23 => Settings::Orthoscopic,
//             24 => Settings::SpecReflect,
//             25 => Settings::SpecPower,
//             26 => Settings::SweepAngle,
//             27 => Settings::SweepSpeed,
//             28 => Settings::DotHydrogens,
//             29 => Settings::DotRadius,
//             30 => Settings::RayTraceFrames,
//             31 => Settings::CacheFrames,
//             32 => Settings::TrimDots,
//             33 => Settings::CullSpheres,
//             34 => Settings::Test1,
//             35 => Settings::Test2,
//             36 => Settings::SurfaceBest,
//             37 => Settings::SurfaceNormal,
//             38 => Settings::SurfaceQuality,
//             39 => Settings::SurfaceProximity,
//             40 => Settings::NormalWorkaround,
//             41 => Settings::StereoAngle,
//             42 => Settings::StereoShift,
//             43 => Settings::LineSmooth,
//             44 => Settings::LineWidth,
//             45 => Settings::HalfBonds,
//             46 => Settings::StickQuality,
//             47 => Settings::StickOverlap,
//             48 => Settings::StickNub,
//             49 => Settings::AllStates,
//             50 => Settings::Pickable,
//             51 => Settings::AutoShowLines,
//             52 => Settings::IdleDelay,
//             53 => Settings::NoIdle,
//             54 => Settings::FastIdle,
//             55 => Settings::SlowIdle,
//             56 => Settings::RockDelay,
//             57 => Settings::DistCounter,
//             58 => Settings::DashLength,
//             59 => Settings::DashGap,
//             60 => Settings::AutoZoom,
//             61 => Settings::Overlay,
//             62 => Settings::Text,
//             63 => Settings::ButtonMode,
//             64 => Settings::Valence,
//             65 => Settings::NonbondedSize,
//             66 => Settings::LabelColor,
//             67 => Settings::RayTraceFog,
//             68 => Settings::SpheroidScale,
//             69 => Settings::RayTraceFogStart,
//             70 => Settings::SpheroidSmooth,
//             71 => Settings::SpheroidFill,
//             72 => Settings::AutoShowNonbonded,
//             73 => Settings::CacheDisplay,
//             74 => Settings::MeshRadius,
//             75 => Settings::BackfaceCull,
//             76 => Settings::Gamma,
//             77 => Settings::DotWidth,
//             78 => Settings::AutoShowSelections,
//             79 => Settings::AutoHideSelections,
//             80 => Settings::SelectionWidth,
//             81 => Settings::SelectionOverlay,
//             82 => Settings::StaticSingletons,
//             84 => Settings::DepthCue,
//             85 => Settings::Specular,
//             86 => Settings::Shininess,
//             87 => Settings::SphereQuality,
//             88 => Settings::Fog,
//             89 => Settings::IsomeshAutoState,
//             90 => Settings::MeshWidth,
//             91 => Settings::CartoonSampling,
//             92 => Settings::CartoonLoopRadius,
//             93 => Settings::CartoonLoopQuality,
//             94 => Settings::CartoonPower,
//             95 => Settings::CartoonPowerB,
//             96 => Settings::CartoonRectLength,
//             97 => Settings::CartoonRectWidth,
//             98 => Settings::InternalGuiWidth,
//             99 => Settings::InternalGui,
//             100 => Settings::CartoonOvalLength,
//             101 => Settings::CartoonOvalWidth,
//             102 => Settings::CartoonOvalQuality,
//             103 => Settings::CartoonTubeRadius,
//             104 => Settings::CartoonTubeQuality,
//             105 => Settings::CartoonDebug,
//             106 => Settings::RibbonWidth,
//             107 => Settings::DashWidth,
//             108 => Settings::DashRadius,
//             109 => Settings::CgoRayWidthScale,
//             110 => Settings::LineRadius,
//             111 => Settings::CartoonRoundHelices,
//             112 => Settings::CartoonRefineNormals,
//             113 => Settings::CartoonFlatSheets,
//             114 => Settings::CartoonSmoothLoops,
//             115 => Settings::CartoonDumbbellLength,
//             116 => Settings::CartoonDumbbellWidth,
//             117 => Settings::CartoonDumbbellRadius,
//             118 => Settings::CartoonFancyHelices,
//             119 => Settings::CartoonFancySheets,
//             120 => Settings::IgnorePdbSegi,
//             121 => Settings::RibbonThrow,
//             122 => Settings::CartoonThrow,
//             123 => Settings::CartoonRefine,
//             124 => Settings::CartoonRefineTips,
//             125 => Settings::CartoonDiscreteColors,
//             126 => Settings::NormalizeCcp4Maps,
//             127 => Settings::SurfacePoor,
//             128 => Settings::InternalFeedback,
//             129 => Settings::CgoLineWidth,
//             130 => Settings::CgoLineRadius,
//             131 => Settings::Logging,
//             132 => Settings::RobustLogs,
//             133 => Settings::LogBoxSelections,
//             134 => Settings::LogConformations,
//             135 => Settings::ValenceSize,
//             136 => Settings::SurfaceMiserable,
//             137 => Settings::RayOpaqueBackground,
//             138 => Settings::Transparency,
//             139 => Settings::RayTexture,
//             140 => Settings::RayTextureSettings,
//             141 => Settings::SuspendUpdates,
//             142 => Settings::FullScreen,
//             143 => Settings::SurfaceMode,
//             144 => Settings::SurfaceColor,
//             145 => Settings::MeshMode,
//             146 => Settings::MeshColor,
//             147 => Settings::AutoIndicateFlags,
//             148 => Settings::SurfaceDebug,
//             149 => Settings::RayImproveShadows,
//             150 => Settings::SmoothColorTriangle,
//             151 => Settings::RayDefaultRenderer,
//             152 => Settings::FieldOfView,
//             153 => Settings::ReflectPower,
//             154 => Settings::PreserveChempyIds,
//             155 => Settings::SphereScale,
//             156 => Settings::TwoSidedLighting,
//             157 => Settings::SecondaryStructure,
//             158 => Settings::AutoRemoveHydrogens,
//             159 => Settings::RaiseExceptions,
//             160 => Settings::StopOnExceptions,
//             161 => Settings::Sculpting,
//             162 => Settings::AutoSculpt,
//             163 => Settings::SculptVdwScale,
//             164 => Settings::SculptVdwScale14,
//             165 => Settings::SculptVdwWeight,
//             166 => Settings::SculptVdwWeight14,
//             167 => Settings::SculptBondWeight,
//             168 => Settings::SculptAnglWeight,
//             169 => Settings::SculptPyraWeight,
//             170 => Settings::SculptPlanWeight,
//             171 => Settings::SculptingCycles,
//             172 => Settings::SphereTransparency,
//             173 => Settings::SphereColor,
//             174 => Settings::SculptFieldMask,
//             175 => Settings::SculptHbOverlap,
//             176 => Settings::SculptHbOverlapBase,
//             177 => Settings::LegacyVdwRadii,
//             178 => Settings::SculptMemory,
//             179 => Settings::ConnectMode,
//             180 => Settings::CartoonCylindricalHelices,
//             181 => Settings::CartoonHelixRadius,
//             182 => Settings::ConnectCutoff,
//             183 => Settings::SavePdbSs,
//             184 => Settings::SculptLineWeight,
//             185 => Settings::FitIterations,
//             186 => Settings::FitTolerance,
//             187 => Settings::BatchPrefix,
//             188 => Settings::StereoMode,
//             189 => Settings::CgoSphereQuality,
//             190 => Settings::PdbLiteralNames,
//             191 => Settings::WrapOutput,
//             192 => Settings::FogStart,
//             193 => Settings::State,
//             194 => Settings::Frame,
//             195 => Settings::RayShadow,
//             196 => Settings::RibbonTraceAtoms,
//             197 => Settings::Security,
//             198 => Settings::StickTransparency,
//             199 => Settings::RayTransparencyShadows,
//             200 => Settings::SessionVersionCheck,
//             201 => Settings::RayTransparencySpecular,
//             202 => Settings::StereDoublePumpMono,
//             203 => Settings::SphereSolvent,
//             204 => Settings::MeshQuality,
//             205 => Settings::MeshSolvent,
//             206 => Settings::DotSolvent,
//             207 => Settings::RayShadowFudge,
//             208 => Settings::RayTriangleFudge,
//             209 => Settings::DebugPick,
//             210 => Settings::DotColor,
//             211 => Settings::MouseLimit,
//             212 => Settings::MouseScale,
//             213 => Settings::TransparencyMode,
//             214 => Settings::ClampColors,
//             215 => Settings::PymolSpaceMaxRed,
//             216 => Settings::PymolSpaceMaxGreen,
//             217 => Settings::PymolSpaceMaxBlue,
//             218 => Settings::PymolSpaceMinFactor,
//             219 => Settings::RovingOrigin,
//             220 => Settings::RovingLines,
//             221 => Settings::RovingSticks,
//             222 => Settings::RovingSpheres,
//             223 => Settings::RovingLabels,
//             224 => Settings::RovingDelay,
//             225 => Settings::RovingSelection,
//             226 => Settings::RovingByres,
//             227 => Settings::RovingRibbon,
//             228 => Settings::RovingCartoon,
//             229 => Settings::RovingPolarContacts,
//             230 => Settings::RovingPolarCutoff,
//             231 => Settings::RovingNonbonded,
//             232 => Settings::FloatLabels,
//             233 => Settings::RovingDetail,
//             234 => Settings::RovingNbSpheres,
//             235 => Settings::RibbonColor,
//             236 => Settings::CartoonColor,
//             237 => Settings::RibbonSmooth,
//             238 => Settings::AutoColor,
//             239 => Settings::AutoColorNext,
//             240 => Settings::RayInteriorColor,
//             241 => Settings::CartoonHighlightColor,
//             242 => Settings::CoulombUnitsFactor,
//             243 => Settings::CoulombDielectric,
//             244 => Settings::RayInteriorShadows,
//             245 => Settings::RayInteriorTexture,
//             246 => Settings::RovingMap1Name,
//             247 => Settings::RovingMap2Name,
//             248 => Settings::RovingMap3Name,
//             249 => Settings::RovingMap1Level,
//             250 => Settings::RovingMap2Level,
//             251 => Settings::RovingMap3Level,
//             252 => Settings::RovingIsomesh,
//             253 => Settings::RovingIsosurface,
//             254 => Settings::ScenesChanged,
//             255 => Settings::GaussianBAdjust,
//             256 => Settings::PdbStandardOrder,
//             257 => Settings::CartoonSmoothFirst,
//             258 => Settings::CartoonSmoothLast,
//             259 => Settings::CartoonSmoothCycles,
//             260 => Settings::CartoonFlatCycles,
//             261 => Settings::MaxThreads,
//             262 => Settings::ShowProgress,
//             263 => Settings::UseDisplayLists,
//             264 => Settings::CacheMemory,
//             265 => Settings::SimplifyDisplayLists,
//             266 => Settings::RetainOrder,
//             267 => Settings::PdbHetatmSort,
//             268 => Settings::PdbUseTerRecords,
//             269 => Settings::CartoonTraceAtoms,
//             270 => Settings::RayOversampleCutoff,
//             271 => Settings::GaussianResolution,
//             272 => Settings::GaussianBFloor,
//             273 => Settings::SculptNbInterval,
//             274 => Settings::SculptTorsWeight,
//             275 => Settings::SculptTorsTolerance,
//             276 => Settings::StickBall,
//             277 => Settings::StickBallRatio,
//             278 => Settings::StickFixedRadius,
//             279 => Settings::CartoonTransparency,
//             280 => Settings::DashRoundEnds,
//             281 => Settings::HBondMaxAngle,
//             282 => Settings::HBondCutoffCenter,
//             283 => Settings::HBondCutoffEdge,
//             284 => Settings::HBondPowerA,
//             285 => Settings::HBondPowerB,
//             286 => Settings::HBondCone,
//             287 => Settings::SsHelixPsiTarget,
//             288 => Settings::SsHelixPsiInclude,
//             289 => Settings::SsHelixPsiExclude,
//             290 => Settings::SsHelixPhiTarget,
//             291 => Settings::SsHelixPhiInclude,
//             292 => Settings::SsHelixPhiExclude,
//             293 => Settings::SsStrandPsiTarget,
//             294 => Settings::SsStrandPsiInclude,
//             295 => Settings::SsStrandPsiExclude,
//             296 => Settings::SsStrandPhiTarget,
//             297 => Settings::SsStrandPhiInclude,
//             298 => Settings::SsStrandPhiExclude,
//             299 => Settings::MovieLoop,
//             300 => Settings::PdbRetainIds,
//             301 => Settings::PdbNoEndRecord,
//             302 => Settings::CgoDotWidth,
//             303 => Settings::CgoDotRadius,
//             304 => Settings::DeferUpdates,
//             305 => Settings::NormalizeOMaps,
//             306 => Settings::SwapDsn6Bytes,
//             307 => Settings::PdbInsertionsGoFirst,
//             308 => Settings::RovingOriginZ,
//             309 => Settings::RovingOriginZCushion,
//             310 => Settings::SpecularIntensity,
//             311 => Settings::OverlayLines,
//             312 => Settings::RayTransparencySpecCut,
//             313 => Settings::InternalPrompt,
//             314 => Settings::Normal,
//             315 => Settings::RayBlendColors,
//             316 => Settings::RayBlendRed,
//             317 => Settings::RayBlendGreen,
//             318 => Settings::RayBlendBlue,
//             319 => Settings::PngScreenGamma,
//             320 => Settings::PngFileGamma,
//             321 => Settings::EditorLabelFragments,
//             322 => Settings::InternalGuiControlSize,
//             323 => Settings::AutoDss,
//             324 => Settings::TransparencyPickingMode,
//             325 => Settings::VirtualTrackball,
//             326 => Settings::PdbReformatNamesMode,
//             327 => Settings::RayPixelScale,
//             328 => Settings::LabelFontId,
//             329 => Settings::PdbConectAll,
//             330 => Settings::ButtonModeName,
//             331 => Settings::SurfaceType,
//             332 => Settings::DotNormals,
//             333 => Settings::SessionMigration,
//             334 => Settings::MeshNormals,
//             335 => Settings::MeshType,
//             336 => Settings::DotLighting,
//             337 => Settings::MeshLighting,
//             338 => Settings::SurfaceSolvent,
//             339 => Settings::TriangleMaxPasses,
//             340 => Settings::RayInteriorReflect,
//             341 => Settings::InternalGuiMode,
//             342 => Settings::SurfaceCarveSelection,
//             343 => Settings::SurfaceCarveState,
//             344 => Settings::SurfaceCarveCutoff,
//             345 => Settings::SurfaceClearSelection,
//             346 => Settings::SurfaceClearState,
//             347 => Settings::SurfaceClearCutoff,
//             348 => Settings::SurfaceTrimCutoff,
//             349 => Settings::SurfaceTrimFactor,
//             350 => Settings::RayMaxPasses,
//             351 => Settings::ActiveSelections,
//             352 => Settings::RayTransparencyContrast,
//             353 => Settings::SeqView,
//             354 => Settings::MouseSelectionMode,
//             355 => Settings::SeqViewLabelSpacing,
//             356 => Settings::SeqViewLabelStart,
//             357 => Settings::SeqViewFormat,
//             358 => Settings::SeqViewLocation,
//             359 => Settings::SeqViewOverlay,
//             360 => Settings::AutoClassifyAtoms,
//             361 => Settings::CartoonNucleicAcidMode,
//             362 => Settings::SeqViewColor,
//             363 => Settings::SeqViewLabelMode,
//             364 => Settings::SurfaceRampAboveMode,
//             365 => Settings::Stereo,
//             366 => Settings::WizardPromptMode,
//             367 => Settings::CoulombCutoff,
//             368 => Settings::SliceTrackCamera,
//             369 => Settings::SliceHeightScale,
//             370 => Settings::SliceHeightMap,
//             371 => Settings::SliceGrid,
//             372 => Settings::SliceDynamicGrid,
//             373 => Settings::SliceDynamicGridResolution,
//             374 => Settings::PdbInsureOrthogonal,
//             375 => Settings::RayDirectShade,
//             376 => Settings::StickColor,
//             377 => Settings::CartoonPuttyRadius,
//             378 => Settings::CartoonPuttyQuality,
//             379 => Settings::CartoonPuttyScaleMin,
//             380 => Settings::CartoonPuttyScaleMax,
//             381 => Settings::CartoonPuttyScalePower,
//             382 => Settings::CartoonPuttyRange,
//             383 => Settings::CartoonSideChainHelper,
//             384 => Settings::SurfaceOptimizeSubsets,
//             385 => Settings::Multiplex,
//             386 => Settings::TextureFonts,
//             387 => Settings::PqrNoChainId,
//             388 => Settings::Animation,
//             389 => Settings::AnimationDuration,
//             390 => Settings::SceneAnimation,
//             391 => Settings::LineStickHelper,
//             392 => Settings::RayOrthoscopic,
//             393 => Settings::RibbonSideChainHelper,
//             394 => Settings::SelectionWidthMax,
//             395 => Settings::SelectionWidthScale,
//             396 => Settings::SceneCurrentName,
//             397 => Settings::Presentation,
//             398 => Settings::PresentationMode,
//             399 => Settings::PdbTruncateResidueName,
//             400 => Settings::SceneLoop,
//             401 => Settings::SweepMode,
//             402 => Settings::SweepPhase,
//             403 => Settings::SceneRestartMovieDelay,
//             404 => Settings::MouseRestartMovieDelay,
//             405 => Settings::AngleSize,
//             406 => Settings::AngleLabelPosition,
//             407 => Settings::DihedralSize,
//             408 => Settings::DihedralLabelPosition,
//             409 => Settings::DeferBuildsMode,
//             410 => Settings::SeqViewDiscreteByState,
//             411 => Settings::SceneAnimationDuration,
//             412 => Settings::Wildcard,
//             413 => Settings::AtomNameWildcard,
//             414 => Settings::IgnoreCase,
//             415 => Settings::PresentationAutoQuit,
//             416 => Settings::EditorAutoDihedral,
//             417 => Settings::PresentationAutoStart,
//             418 => Settings::ValidateObjectNames,
//             419 => Settings::UnusedBooleanDefTrue,
//             420 => Settings::AutoShowSpheres,
//             421 => Settings::SphereMode,
//             422 => Settings::SpherePointMaxSize,
//             423 => Settings::SpherePointSize,
//             424 => Settings::PdbHonorModelNumber,
//             425 => Settings::RankAssistedSorts,
//             426 => Settings::RibbonNucleicAcidMode,
//             427 => Settings::CartoonRingMode,
//             428 => Settings::CartoonRingWidth,
//             429 => Settings::CartoonRingColor,
//             430 => Settings::CartoonRingFinder,
//             431 => Settings::CartoonTubeCap,
//             432 => Settings::CartoonLoopCap,
//             433 => Settings::NvidiaBugs,
//             434 => Settings::ImageDotsPerInch,
//             435 => Settings::OpaqueBackground,
//             436 => Settings::DrawFrames,
//             437 => Settings::ShowAlphaChecker,
//             438 => Settings::MatrixMode,
//             439 => Settings::EditorAutoOrigin,
//             440 => Settings::SessionFile,
//             441 => Settings::CgoTransparency,
//             442 => Settings::LegacyMouseZoom,
//             443 => Settings::AutoNumberSelections,
//             444 => Settings::SculptVdwVisMode,
//             445 => Settings::SculptVdwVisMin,
//             446 => Settings::SculptVdwVisMid,
//             447 => Settings::SculptVdwVisMax,
//             448 => Settings::CartoonLadderMode,
//             449 => Settings::CartoonLadderRadius,
//             450 => Settings::CartoonLadderColor,
//             451 => Settings::CartoonNucleicAcidColor,
//             452 => Settings::CartoonRingTransparency,
//             453 => Settings::LabelSize,
//             454 => Settings::SpecDirect,
//             455 => Settings::LightCount,
//             456 => Settings::Light2,
//             457 => Settings::Light3,
//             458 => Settings::HideUnderscoreNames,
//             459 => Settings::SelectionRoundPoints,
//             460 => Settings::DistanceExclusion,
//             461 => Settings::HBondExclusion,
//             462 => Settings::LabelShadowMode,
//             463 => Settings::Light4,
//             464 => Settings::Light5,
//             465 => Settings::Light6,
//             466 => Settings::Light7,
//             467 => Settings::LabelOutlineColor,
//             468 => Settings::RayTraceMode,
//             469 => Settings::RayTraceGain,
//             470 => Settings::SelectionVisibleOnly,
//             471 => Settings::LabelPosition,
//             472 => Settings::RayTraceDepthFactor,
//             473 => Settings::RayTraceSlopeFactor,
//             474 => Settings::RayTraceDiscoFactor,
//             475 => Settings::RayShadowDecayFactor,
//             476 => Settings::RayInteriorMode,
//             477 => Settings::RayLegacyLighting,
//             478 => Settings::SculptAutoCenter,
//             479 => Settings::PdbDiscreteChains,
//             480 => Settings::PdbUnbondCations,
//             481 => Settings::SculptTriScale,
//             482 => Settings::SculptTriWeight,
//             483 => Settings::SculptTriMin,
//             484 => Settings::SculptTriMax,
//             485 => Settings::SculptTriMode,
//             486 => Settings::PdbEchoTags,
//             487 => Settings::ConnectBonded,
//             488 => Settings::SpecDirectPower,
//             489 => Settings::Light8,
//             490 => Settings::Light9,
//             491 => Settings::RayShadowDecayRange,
//             492 => Settings::SpecCount,
//             493 => Settings::SculptMinScale,
//             494 => Settings::SculptMinWeight,
//             495 => Settings::SculptMinMin,
//             496 => Settings::SculptMinMax,
//             497 => Settings::SculptMaxScale,
//             498 => Settings::SculptMaxWeight,
//             499 => Settings::SculptMaxMin,
//             500 => Settings::SculptMaxMax,
//             501 => Settings::SurfaceCircumscribe,
//             502 => Settings::SculptAvdWeight,
//             503 => Settings::SculptAvdGap,
//             504 => Settings::SculptAvdRange,
//             505 => Settings::SculptAvdExcl,
//             506 => Settings::AsyncBuilds,
//             507 => Settings::FetchPath,
//             508 => Settings::CartoonRingRadius,
//             509 => Settings::RayColorRamps,
//             510 => Settings::RayHintCamera,
//             511 => Settings::RayHintShadow,
//             512 => Settings::StickValenceScale,
//             513 => Settings::SeqViewAlignment,
//             514 => Settings::SeqViewUnalignedMode,
//             515 => Settings::SeqViewUnalignedColor,
//             516 => Settings::SeqViewFillChar,
//             517 => Settings::SeqViewFillColor,
//             518 => Settings::SeqViewLabelColor,
//             519 => Settings::SurfaceCarveNormalCutoff,
//             520 => Settings::TraceAtomsMode,
//             521 => Settings::SessionChanged,
//             522 => Settings::RayClipShadows,
//             523 => Settings::MouseWheelScale,
//             524 => Settings::NonbondedTransparency,
//             525 => Settings::RaySpecLocal,
//             526 => Settings::LineColor,
//             527 => Settings::RayLabelSpecular,
//             528 => Settings::MeshSkip,
//             529 => Settings::LabelDigits,
//             530 => Settings::LabelDistanceDigits,
//             531 => Settings::LabelAngleDigits,
//             532 => Settings::LabelDihedralDigits,
//             533 => Settings::SurfaceNegativeVisible,
//             534 => Settings::SurfaceNegativeColor,
//             535 => Settings::MeshNegativeVisible,
//             536 => Settings::MeshNegativeColor,
//             537 => Settings::GroupAutoMode,
//             538 => Settings::GroupFullMemberNames,
//             539 => Settings::GradientMaxLength,
//             540 => Settings::GradientMinLength,
//             541 => Settings::GradientMinSlope,
//             542 => Settings::GradientNormalMinDot,
//             543 => Settings::GradientStepSize,
//             544 => Settings::GradientSpacing,
//             545 => Settings::GradientSymmetry,
//             546 => Settings::RayTraceColor,
//             547 => Settings::GroupArrowPrefix,
//             548 => Settings::SuppressHidden,
//             549 => Settings::SessionCompression,
//             550 => Settings::MovieFps,
//             551 => Settings::RayTransparencyOblique,
//             552 => Settings::RayTraceTransCutoff,
//             553 => Settings::RayTracePersistCutoff,
//             554 => Settings::RayTransparencyObliquePower,
//             555 => Settings::RayScatter,
//             556 => Settings::HBondFromProton,
//             557 => Settings::AutoCopyImages,
//             558 => Settings::MoeSeparateChains,
//             559 => Settings::TransparencyGlobalSort,
//             560 => Settings::HideLongBonds,
//             561 => Settings::AutoRenameDuplicateObjects,
//             562 => Settings::PdbHetatmGuessValences,
//             563 => Settings::EllipsoidQuality,
//             564 => Settings::CgoEllipsoidQuality,
//             565 => Settings::MovieAnimateByFrame,
//             566 => Settings::RampBlendNearbyColors,
//             567 => Settings::AutoDeferBuilds,
//             568 => Settings::EllipsoidProbability,
//             569 => Settings::EllipsoidScale,
//             570 => Settings::EllipsoidColor,
//             571 => Settings::EllipsoidTransparency,
//             572 => Settings::MovieRock,
//             573 => Settings::CacheMode,
//             574 => Settings::DashColor,
//             575 => Settings::AngleColor,
//             576 => Settings::DihedralColor,
//             577 => Settings::GridMode,
//             578 => Settings::CacheMax,
//             579 => Settings::GridSlot,
//             580 => Settings::GridMax,
//             581 => Settings::CartoonPuttyTransform,
//             582 => Settings::Rock,
//             583 => Settings::ConeQuality,
//             584 => Settings::PdbFormalCharges,
//             585 => Settings::AtiBugs,
//             586 => Settings::GeometryExportMode,
//             587 => Settings::MouseGrid,
//             588 => Settings::MeshCutoff,
//             589 => Settings::MeshCarveSelection,
//             590 => Settings::MeshCarveState,
//             591 => Settings::MeshCarveCutoff,
//             592 => Settings::MeshClearSelection,
//             593 => Settings::MeshClearState,
//             594 => Settings::MeshClearCutoff,
//             595 => Settings::MeshGridMax,
//             596 => Settings::SessionCacheOptimize,
//             597 => Settings::SdofDragScale,
//             598 => Settings::SceneButtonsMode,
//             599 => Settings::SceneButtons,
//             600 => Settings::MapAutoExpandSym,
//             601 => Settings::ImageCopyAlways,
//             602 => Settings::MaxUps,
//             603 => Settings::AutoOverlay,
//             604 => Settings::StickBallColor,
//             605 => Settings::StickHScale,
//             606 => Settings::SculptPyraInvWeight,
//             607 => Settings::KeepAlive,
//             608 => Settings::FitKabsch,
//             609 => Settings::StereoDynamicStrength,
//             610 => Settings::DynamicWidth,
//             611 => Settings::DynamicWidthFactor,
//             612 => Settings::DynamicWidthMin,
//             613 => Settings::DynamicWidthMax,
//             614 => Settings::DrawMode,
//             615 => Settings::CleanElectroMode,
//             616 => Settings::ValenceMode,
//             617 => Settings::ShowFrameRate,
//             618 => Settings::MoviePanel,
//             619 => Settings::MouseZScale,
//             620 => Settings::MovieAutoStore,
//             621 => Settings::MovieAutoInterpolate,
//             622 => Settings::MoviePanelRowHeight,
//             623 => Settings::SceneFrameMode,
//             624 => Settings::SurfaceCavityMode,
//             625 => Settings::SurfaceCavityRadius,
//             626 => Settings::SurfaceCavityCutoff,
//             627 => Settings::MotionPower,
//             628 => Settings::MotionBias,
//             629 => Settings::MotionSimple,
//             630 => Settings::MotionLinear,
//             631 => Settings::MotionHand,
//             632 => Settings::PdbIgnoreConect,
//             633 => Settings::EditorBondCycleMode,
//             634 => Settings::MovieQuality,
//             635 => Settings::LabelAnchor,
//             636 => Settings::FetchHost,
//             637 => Settings::DynamicMeasures,
//             638 => Settings::NeighborCutoff,
//             639 => Settings::HeavyNeighborCutoff,
//             640 => Settings::PolarNeighborCutoff,
//             641 => Settings::SurfaceResidueCutoff,
//             642 => Settings::SurfaceUseShader,
//             643 => Settings::CartoonUseShader,
//             644 => Settings::StickUseShader,
//             645 => Settings::LineUseShader,
//             646 => Settings::SphereUseShader,
//             647 => Settings::UseShaders,
//             648 => Settings::ShadersFromDisk,
//             649 => Settings::VolumeBitDepth,
//             650 => Settings::VolumeColor,
//             651 => Settings::VolumeLayers,
//             652 => Settings::VolumeDataRange,
//             653 => Settings::AutoDeferAtomCount,
//             654 => Settings::DefaultRefmacNames,
//             655 => Settings::DefaultPhenixNames,
//             656 => Settings::DefaultPhenixNoFillNames,
//             657 => Settings::DefaultBusterNames,
//             658 => Settings::DefaultFofcMapRep,
//             659 => Settings::Default2fofcMapRep,
//             660 => Settings::AtomTypeFormat,
//             661 => Settings::AutocloseDialogs,
//             662 => Settings::BgGradient,
//             663 => Settings::BgRgbTop,
//             664 => Settings::BgRgbBottom,
//             665 => Settings::RayVolume,
//             666 => Settings::RibbonTransparency,
//             667 => Settings::StateCounterMode,
//             668 => Settings::CgoUseShader,
//             669 => Settings::CgoShaderUbColor,
//             670 => Settings::CgoShaderUbNormal,
//             671 => Settings::CgoLighting,
//             672 => Settings::MeshUseShader,
//             673 => Settings::StickDebug,
//             674 => Settings::CgoDebug,
//             675 => Settings::StickRoundNub,
//             676 => Settings::StickGoodGeometry,
//             677 => Settings::StickAsCylinders,
//             678 => Settings::MeshAsCylinders,
//             679 => Settings::LineAsCylinders,
//             680 => Settings::RibbonAsCylinders,
//             681 => Settings::RibbonUseShader,
//             682 => Settings::ExclDisplayListsShaders,
//             683 => Settings::DashUseShader,
//             684 => Settings::DashAsCylinders,
//             685 => Settings::NonbondedUseShader,
//             686 => Settings::NonbondedAsCylinders,
//             687 => Settings::CylindersShaderFilterFaces,
//             688 => Settings::NbSpheresSize,
//             689 => Settings::NbSpheresQuality,
//             690 => Settings::NbSpheresUseShader,
//             691 => Settings::RenderAsCylinders,
//             692 => Settings::AlignmentAsCylinders,
//             693 => Settings::CartoonNucleicAcidAsCylinders,
//             694 => Settings::CgoShaderUbFlags,
//             695 => Settings::AntialiasShader,
//             696 => Settings::OffscreenRenderingMultiplier,
//             697 => Settings::CylinderShaderFfWorkaround,
//             698 => Settings::SurfaceColorSmoothing,
//             699 => Settings::SurfaceColorSmoothingThreshold,
//             700 => Settings::DotUseShader,
//             701 => Settings::DotAsSheres,
//             702 => Settings::AmbientOcclusionMode,
//             703 => Settings::AmbientOcclusionScale,
//             704 => Settings::AmbientOcclusionSmooth,
//             705 => Settings::SmoothHalfBonds,
//             706 => Settings::AnaglyphMode,
//             707 => Settings::EditLight,
//             708 => Settings::SuspendUndo,
//             709 => Settings::SuspendUndoAtomCount,
//             710 => Settings::SuspendDeferred,
//             711 => Settings::PickSurface,
//             712 => Settings::BgImageFilename,
//             713 => Settings::BgImageMode,
//             714 => Settings::BgImageTilesize,
//             715 => Settings::BgImageLinear,
//             716 => Settings::LoadObjectPropsDefault,
//             717 => Settings::LoadAtomPropsDefault,
//             718 => Settings::LabelPlacementOffset,
//             719 => Settings::PdbConectNodup,
//             720 => Settings::LabelConnector,
//             721 => Settings::LabelConnectorMode,
//             722 => Settings::LabelConnectorColor,
//             723 => Settings::LabelConnectorWidth,
//             724 => Settings::LabelConnectorExtLength,
//             725 => Settings::LabelBgColor,
//             726 => Settings::UseGeometryShaders,
//             727 => Settings::LabelRelativeMode,
//             728 => Settings::LabelScreenPoint,
//             729 => Settings::LabelMultilineSpacing,
//             730 => Settings::LabelMultilineJustification,
//             731 => Settings::LabelPadding,
//             732 => Settings::LabelBgTransparency,
//             733 => Settings::LabelBgOutline,
//             734 => Settings::RayLabelConnectorFlat,
//             735 => Settings::DashTransparency,
//             736 => Settings::PickLabels,
//             737 => Settings::LabelZTarget,
//             738 => Settings::SessionEmbedsData,
//             739 => Settings::VolumeMode,
//             740 => Settings::Trilines,
//             741 => Settings::ColladaExportLighting,
//             742 => Settings::ColladaGeometryMode,
//             743 => Settings::PrecomputedLighting,
//             744 => Settings::Chromadepth,
//             745 => Settings::PseExportVersion,
//             746 => Settings::CifUseAuth,
//             747 => Settings::Assembly,
//             748 => Settings::CifKeepinmemory,
//             749 => Settings::PseBinaryDump,
//             750 => Settings::CartoonGapCutoff,
//             751 => Settings::IgnoreCaseChain,
//             752 => Settings::ValenceZeroScale,
//             753 => Settings::ValenceZeroMode,
//             754 => Settings::AutoShowClassified,
//             755 => Settings::ColladaBackgroundBox,
//             756 => Settings::Pick32bit,
//             757 => Settings::CartoonAllAlt,
//             758 => Settings::DisplayScaleFactor,
//             759 => Settings::PickShading,
//             760 => Settings::FetchTypeDefault,
//             761 => Settings::EditorAutoMeasure,
//             762 => Settings::SurfaceSmoothEdges,
//             763 => Settings::ChemCompCartnUse,
//             764 => Settings::ColoredFeedback,
//             765 => Settings::SdfWriteZeroOrderBonds,
//             766 => Settings::CifMetalcAsZeroOrderBonds,
//             767 => Settings::SeqViewGapMode,
//             768 => Settings::InternalGuiNameColorMode,
//             769 => Settings::OpenvrGuiFov,
//             770 => Settings::OpenvrGuiAlpha,
//             771 => Settings::OpenvrGuiUseAlpha,
//             772 => Settings::OpenvrGuiSceneColor,
//             773 => Settings::OpenvrGuiSceneAlpha,
//             774 => Settings::OpenvrGuiBackColor,
//             775 => Settings::OpenvrGuiBackAlpha,
//             776 => Settings::OpenvrGuiUseBackdrop,
//             777 => Settings::OpenvrGuiOverlay,
//             778 => Settings::OpenvrGuiText,
//             779 => Settings::OpenvrDisableClipping,
//             780 => Settings::OpenvrNearPlane,
//             781 => Settings::OpenvrFarPlane,
//             782 => Settings::OpenvrCutLaser,
//             783 => Settings::OpenvrLaserWidth,
//             784 => Settings::OpenvrGuiDistance,
//             785 => Settings::CartoonSmoothCylinderCycles,
//             786 => Settings::CartoonSmoothCylinderWindow,
//             787 => Settings::IsosurfaceAlgorithm,
//             788 => Settings::CellCentered,
//             789 => Settings::HalogenBondDistance,
//             790 => Settings::HalogenBondAsDonorMinDonorAngle,
//             791 => Settings::HalogenBondAsDonorMinAcceptorAngle,
//             792 => Settings::HalogenBondAsAcceptorMinDonorAngle,
//             793 => Settings::HalogenBondAsAcceptorMinAcceptorAngle,
//             794 => Settings::HalogenBondAsAcceptorMaxAcceptorAngle,
//             795 => Settings::SaltBridgeDistance,
//             796 => Settings::UseTessellationShaders,
//             797 => Settings::CellColor,
//             _ => panic!("Invalid Settings value: {}", value),
//         }
//     }
// }

// Implement Into<i32> for Settings
// impl Into<i32> for Settings {
//     fn into(self) -> i32 {
//         self as i32
//     }
// }

#[derive(Debug, Deserialize, Serialize)]
// #[serde(from = "(SettingsEnum, i32, CustomValue)")]
pub struct Settings {
    setting: SettingsEnum,
    label: i32,
    value: CustomValue,
}
// // Implement the conversion from a tuple to your Settings struct
// impl From<(SettingsEnum, i32, CustomValue)> for Settings {
//     fn from(tuple: (SettingsEnum, i32, CustomValue)) -> Self {
//         Settings {
//             // setting: SettingsEnum::from(tuple.0),
//             setting: tuple.0,
//             label: tuple.1,
//             value: tuple.2,
//         }
//     }
// }

// Todo:
//
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
