use crate::molviewspec::nodes::{ComponentExpression, ComponentSelector};
use itertools::Itertools;
use pdbtbx::{self, Residue, PDB};
use serde::{Deserialize, Deserializer, Serialize};
use serde_pickle::{from_value, Value};

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

/// PyObjectMolecule:
///
/// - [pymol code](https://github.com/schrodinger/pymol-open-source/blob/master/layer2/ObjectMolecule2.cpp#L3524)
/// - [ObjectMolecule](https://github.com/schrodinger/pymol-open-source/blob/03d7a7fcf0bd95cd93d710a1268dbace2ed77765/layer2/ObjectMolecule.h#L58)
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
    // /// https://github.com/schrodinger/pymol-open-source/blob/master/layer2/ObjectMolecule2.cpp#L3037
    pub bond: Vec<Bond>,
    // /// https://github.com/schrodinger/pymol-open-source/blob/master/layer2/ObjectMolecule2.cpp#L3248
    // /// https://github.com/schrodinger/pymol-open-source/blob/master/layer2/AtomInfo.cpp#L792
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
    pub fn get_str(&self) -> String {
        // get_str
        // https://github.com/schrodinger/pymol-open-source/blob/03d7a7fcf0bd95cd93d710a1268dbace2ed77765/layer4/Cmd.cpp#L3877
        // PymolMoleculeExporter
        // https://github.com/schrodinger/pymol-open-source/blob/master/layer3/MoleculeExporter.cpp#L1627
        // PDB Exporter
        // MoleculeExporterPDB
        // https://github.com/schrodinger/pymol-open-source/blob/master/layer3/MoleculeExporter.cpp#L439
        //
        //
        // exporter->init(G);
        //
        // exporter->setMulti(multi);
        // exporter->setRefObject(ref_object, ref_state);
        // exporter->execute(sele, state);
        //
        // CoordSetAtomToPDBStrVLA
        //
        // Variables:
        //  - m_iter: SeleCoordIterator m_iter;
        //
        // // Atoms:
        //     // how do we check for multiple objects?
        //     // m_iter.obj defines the object. By number? by name?
        //     - iterate through the coordinates
        //     - check for multi

        // update transformation matrices
        // updateMatrix(m_mat_full, true);
        // updateMatrix(m_mat_move, false);

        // beginCoordSet();
        // m_last_cs = m_iter.cs;

        // for bonds
        // if (!m_tmpids[m_iter.getAtm()]) {
        //   m_id = m_retain_ids ? m_iter.getAtomInfo()->id : (m_id + 1);
        //   m_tmpids[m_iter.getAtm()] = m_id;
        // }

        //

        // for coord in &self.coord_set {
        //     println!("{:?}", coord);
        // }
        "test".to_string()
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
    ///
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
        // println!("{:?}", model);

        // Create PDB from Models
        let mut pdb = PDB::new();
        pdb.add_model(model);

        // Add Bonds Here
        for bond in &self.bond {
            // println!("{:?}", bond);
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

#[derive(Debug, Deserialize, Serialize)]
pub struct SessionSelectorList(Vec<SessionSelector>);
impl SessionSelectorList {
    pub fn get_selectors(&self) -> &Vec<SessionSelector> {
        &self.0
    }
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
                label_entity_id: Some(self.id.clone()),
                atom_id: Some(idx32),
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

/// Coord Set
///
/// - [pymol_coordset](https://github.com/schrodinger/pymol-open-source/blob/master/layer2/CoordSet.cpp#L363)
/// - [pymol_coordset_settings]( https://github.com/schrodinger/pymol-open-source/blob/master/layer1/Setting.cpp#L962)
#[derive(Debug, Deserialize, Serialize)]
pub struct CoordSet {
    pub n_index: i32,         // 1519
    n_at_index: i32,          // 1519
    pub coord: Vec<f32>,      // len -== 4556 ( 1519 *3 )
    pub idx_to_atm: Vec<i32>, // 1 - 1518
    pub atm_to_idx: Option<Vec<i32>>,
    pub name: String,
    pub setting: Vec<Option<bool>>, // punting on this
    pub lab_pos: Option<bool>,      // might be wrong
    field_9: Option<bool>,          // probably wrong...
    // /// https://github.com/schrodinger/pymol-open-source/blob/master/layer1/CGO.cpp#L220
    pub sculpt_cgo: Option<(i32, Vec<f32>)>,           //
    pub atom_state_settings: Option<Vec<Option<i32>>>, //
    /// [symettry_settings](https://github.com/schrodinger/pymol-open-source/blob/master/layer1/Symmetry.cpp#L30)
    pub symmetry: Option<Vec<(((i32, i32, i32), (i32, i32, i32)), String)>>,
}

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

#[derive(Debug, Deserialize, Serialize)]
pub struct Bond {
    pub index_1: i32,
    pub index_2: i32,
    pub order: i32,
    pub id: i32,
    pub stereo: i32,
    pub unique_id: i32,
    pub has_setting: i32,
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
