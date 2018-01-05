from __future__ import print_function

__title__ = "FreeCAD ANSYS rst reader based on pyansys"
__author__ = "AMStuff"
__url__ = "http://www.freecadweb.org"

## @package importAnsysRstMesh
#  \ingroup FEM
#  \brief FreeCAD Ansys rst result reader for FEM workbench

import FreeCAD
import Fem
import os

########## generic FreeCAD import and export methods ##########
if open.__module__ == '__builtin__':
    # because we'll redefine open below (Python2)
    pyopen = open
elif open.__module__ == 'io':
    # because we'll redefine open below (Python3)
    pyopen = open


def open(filename):
    """
    called when freecad opens a file
    """
    docname = os.path.splitext(os.path.basename(filename))[0]
    insert(filename, docname)


def insert(filename, docname):
    """
    called when freecad wants to import a file
    """
    try:
        doc = FreeCAD.getDocument(docname)
    except NameError:
        doc = FreeCAD.newDocument(docname)
    FreeCAD.ActiveDocument = doc
    import_ansys_rst(filename)


########## module specific methods ##########
def import_ansys_rst(filename):
    """
    insert a FreeCAD FEM Mesh object in the ActiveDocument
    """

    try:
        import pyansys
    except ModuleNotFoundError:
        FreeCAD.Console.PrintError("ANSYS read module requires pyansys to be"
                                   "installed.\n")
        return

    import importToolsFem
    import ObjectsFem
    import tempfile

    pya = pyansys.ResultReader(filename)

    # Make Result Object
    results_name = os.path.basename(os.path.splitext(filename)[0])
    result_obj = ObjectsFem.makeResultMechanical(FreeCAD.ActiveDocument,
                                                 results_name)

    # get geometry by writing a vtu file with pyansys
    tf = tempfile.mkstemp(suffix='.vtu')[1]
    pya.grid.Write(tf)
    # readResult always creates a new femmesh named ResultMesh
    Fem.readResult(tf, result_obj.Name)
    os.remove(tf)

    # calculate the span (thanks to importCcxFrdResults.py)
    nodes = []
    for k, v in result_obj.Mesh.FemMesh.Nodes.items():
        nodes.append(list(v))
    p_x_max, p_y_max, p_z_max = map(max, zip(*nodes))
    p_x_min, p_y_min, p_z_min = map(min, zip(*nodes))
    x_span = abs(p_x_max - p_x_min)
    y_span = abs(p_y_max - p_y_min)
    z_span = abs(p_z_max - p_z_min)
    span = max(x_span, y_span, z_span)

    # get Results
    result_set = get_results(pya)
    results = importToolsFem.fill_femresult_mechanical(result_obj, result_set,
                                                       span)

    # Add extra information
    for key, res in result_set['stress_i'].items():
        attrname = 'NodalStress{:03d}'.format(int(key))
        results.addProperty("App::PropertyVectorList", attrname,
                            "Fem userdefined", attrname, True)
        current_res = [FreeCAD.Vector(b[0:3]) for a in sorted(res.items())
                       for b in list(a)[1:]]
        setattr(results, attrname, current_res)

    for key, res in result_set['nodal_i'].items():
        attrname = 'NodalResult{:03d}'.format(int(key))
        results.addProperty("App::PropertyVectorList", attrname,
                            "Fem userdefined", attrname, True)
        current_res = [FreeCAD.Vector(b[0:3]) for a in sorted(res.items())
                       for b in list(a)[1:]]
        setattr(results, attrname, current_res)


def get_results(pya):
    """
    converts pyansys results into a FreeCAD-compatible format
    """
    result_set = {}
    result_set['stress_i'] = {}
    result_set['nodal_i'] = {}
    # number sets
    if 'nsets' in pya.resultheader:
        nnum = pya.nnum.tolist()
        number_results = pya.resultheader['nsets']

        # First we're retrieving and safing all results
        stress_i = {}
        nodal_i = {}
        for i in range(0, number_results):
            nodalresults = []
            nodalstress = []

            # There's got to be a better way
            try:
                # This will usually be the displacement, but is user defined
                # FIXME is there a way to tell what it is?
                nodalresults = pya.GetNodalResult(i)
            except:
                pass

            try:
                nodalstress = pya.NodalStress(i)
            except:
                pass

            # Sx Sy Sz Sxy Syz Sxz
            if len(nodalstress):
                stress = {}
                for j, vec in enumerate(nodalstress):
                    stress[nnum[j]] = vec.tolist()

                nodalstress = []
                stress_i[i] = stress

            # Nodal Results(can be anything basically)
            if len(nodalresults):
                nodal = {}

                # Convert to FreeCAD Vector if applicable, otherwise convert to
                # list
                if len(nodalresults[0]) == 3:
                    for j, vec in enumerate(nodalresults):
                        nodal[nnum[j]] = FreeCAD.Vector(vec.tolist())
                else:
                    for j, vec in enumerate(nodalresults):
                        nodal[nnum[j]] = vec.tolist()

                nodalresults = []
                nodal_i[i] = nodal

        result_set['stress_i'] = stress_i
        result_set['nodal_i'] = nodal_i

        # convert first applicable stress set to FreeCAD representation
        wrote_stress = False
        wrote_nodal = False
        if stress_i:
            for key, stress_set in stress_i.items():
                # check if there are enough values
                if len(stress_set) == len(nnum):
                    first_dict_entry = next(iter(stress_set))
                    if len(stress_set[first_dict_entry]) == 6:
                        result_set['stress'] = stress_set
                        # Also fill stressv
                        stressv = {}
                        for d, vec in stress_set.items():
                            stressv[d] = FreeCAD.Vector(vec[0:3])
                        result_set['stressv'] = stressv
                        wrote_stress = True
                        break

        if nodal_i:
            for key, nodal_set in nodal_i.items():
                if len(nodal_set) == len(nnum):
                    first_dict_entry = next(iter(nodal_set))
                    if isinstance(nodal_set[first_dict_entry], FreeCAD.Vector):
                        result_set['disp'] = nodal_set
                        wrote_nodal = True
                        break

        # if we have displacements without stresses, fake having stresses to
        # keep freecad happy
        # if wrote_nodal and not wrote_stress:
            # stress = {}
            # stressv = {}
            # for el in nnum:
                # stress[el] = [0,0,0,0,0,0]
                # stressv[el] = FreeCAD.Vector(0,0,0)
            # result_set['stress'] = stress
            # result_set['stressv'] = stressv

    return result_set
