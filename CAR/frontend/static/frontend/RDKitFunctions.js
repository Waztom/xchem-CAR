var Module = {
  onRuntimeInitialized: initRDKitModule().then(function (instance) {
    Module = instance;
    console.log("version: " + Module.version());
  }),
};

function draw(mol) {
  var svg = mol.get_svg();
  if (svg == "") return;
  var ob = document.getElementById("drawing");
  ob.outerHTML = "<div id='drawing'>" + svg + "</div>";
}

function draw_with_highlights(mol, atomids) {
  var svg = mol.get_svg_with_highlights(atomids);
  if (svg == "") return;
  var ob = document.getElementById("drawing");
  ob.outerHTML = "<div id='drawing'>" + svg + "</div>";
}

function callback(text) {
  var mol = Module.get_mol(text);
  if (mol.is_valid()) {
    console.log(mol);
    // draw(mol);
    // var ob = document.getElementById("can_smiles");
    // ob.outerHTML = "<div id='can_smiles'>" + mol.get_smiles() + "</div>";
    // var descrs = JSON.parse(mol.get_descriptors());
    // var db = document.getElementById("descrs");
    // db.outerHTML =
    //   "<div id='descrs'>" +
    //   "<b>AMW:</b> " +
    //   descrs.amw +
    //   "<br /><b>MolLogP:</b> " +
    //   descrs.CrippenClogP +
    //   "<br /><b>MFP2:</b> " +
    //   mol.get_morgan_fp(2, 128) +
    //   "</div>";
  }
  mol.delete();
}

function sma_callback(text) {
  var qmol = Module.get_qmol(text);
  var mol = Module.get_mol(document.getElementById("smiles_input").value);
  if (text == "" && mol.is_valid()) {
    // SMARTS box empty, just draw without highlights
    draw(mol);
  } else if (mol.is_valid() && qmol.is_valid()) {
    var mdetails = mol.get_substruct_match(qmol);
    var match = JSON.parse(mdetails);
    if (match.atoms && match.atoms.length) draw_with_highlights(mol, mdetails);
  }
  mol.delete();
  qmol.delete();
}
