//! PSEData is a struct for loading and serializing pymol PSE data.
//!
//! Currently the parsers are working for small test cases of molecules and selections. Additional parser structs would be required for
//! other PSE data types which include the folloing:
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
use crate::molviewspec::nodes::{self, State};
use crate::pymolparsing::parsing::{
    CustomValue, PyObjectMolecule, PymolSessionObjectData, SessionName,
};
use pdbtbx::PDB;
use serde::{Deserialize, Serialize};
use serde_pickle::de::{from_reader, DeOptions};
use std::io::Read;
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

        let pse_data_vals: serde_pickle::Value = from_reader(&buffer[..], options).unwrap();
        let pse_data: PSEData = serde_pickle::from_value(pse_data_vals).unwrap();
        Ok(pse_data)
    }

    pub fn to_json(&self, file_path: &str) -> Result<(), Box<dyn std::error::Error>> {
        let json = serde_json::to_string_pretty(self)?;
        std::fs::write(file_path, json)?;
        Ok(())
    }

    /// session is where all the action happens
    pub fn get_session_names(&self) -> Vec<String> {
        self.names
            .iter()
            .filter_map(|session_name| {
                session_name
                    .as_ref()
                    .map(|session| session.name.to_string())
            })
            .collect()
    }

    pub fn get_molecule_data(&self) -> Vec<&PyObjectMolecule> {
        self.names
            .iter()
            .filter_map(|session_name| session_name.as_ref())
            .filter_map(|session| match &session.data {
                PymolSessionObjectData::PyObjectMolecule(a) => Some(a),
                _ => None,
            })
            .collect()
    }

    pub fn create_pdb(&self) -> PDB {
        // todo: extend this to more than one molecuelo and/or to modify the global scene
        let moldata = &self.get_molecule_data();
        let first_mol = moldata[0];
        first_mol.to_pdb()
    }

    pub fn save_pdbs(&self, file_path: &str) -> std::io::Result<()> {
        let path = std::path::Path::new(file_path);
        let pdb_folder = path.join("pdb");
        std::fs::create_dir_all(&pdb_folder)?;

        let mut file_list = Vec::new();

        for (index, molecule) in self.get_molecule_data().iter().enumerate() {
            let pdb = molecule.to_pdb();
            let filename = format!("{}.pdb", molecule.get_name());
            let file_path = pdb_folder.join(&filename);

            let _ = pdbtbx::save_pdb(
                &pdb,
                file_path.to_str().expect("Invalid UTF-8 in file path"),
                pdbtbx::StrictnessLevel::Strict,
            );
            file_list.push(filename);
        }

        let contents = file_list.join("\n");
        std::fs::write(path.join("pdb_contents.txt"), contents)?;

        Ok(())
    }

    pub fn create_molviewspec(&self) -> State {
        State::new()
    }
}
