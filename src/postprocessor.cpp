#include "postprocessor.h"
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkTetra.h>
#include <vtkCellArray.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkXMLWriter.h>
#include <vtkSmartPointerBase.h>


postprocessor::postprocessor()
{
    //ctor
}

postprocessor::~postprocessor()
{
    //dtor
}

void postprocessor::output_vtk_mesh(std::string output_location, global_variables &globals,
        domain_geometry &geometry){


    // Create filename.
  stringstream output_filename;
  output_filename <<  output_location + "/grid.vtk";
    ofstream output_file;
    // Open file.
output_file.open(output_filename.str().c_str());

// Write VTK header.
output_file << "# vtk DataFile Version 3.0\n";
  output_file << "fluid_state\n";
output_file << "ASCII\n";
output_file << "DATASET UNSTRUCTURED_GRID\n";
output_file << "DIMENSIONS " << globals. << " " << Ny - 2 << " 1" << "\n";
  output_file << "X_COORDINATES " << Nx << " float\n";
  for(int X = 0; X < Nx; ++X) {
    output_file << X + 0.5 << " ";
}

}
