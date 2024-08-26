#ifndef VTK2DV_H
#define VTK2DV_H

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkWindowToImageFilter.h>
#include <vtkJPEGWriter.h>
#include <vtkPNGWriter.h>
#include <vtkProperty.h>
#include <cstring>
#include <vtkRenderWindowInteractor.h>
#include <vtkLookupTable.h>
#include <vtkColorTransferFunction.h>
#include <vtkLight.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkConeSource.h>
#include <vtkCellData.h>
#include <vtkTimerLog.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>


#include "vecInoue.h"
#include "class2dv.h"
#include "Organ2dv.h"
#include "Cell2dv.h"
#include "Line2dv.h"
#include "Vertex2dv.h"

namespace autoVTK{
    void VTKLine(std::string);
    void VTKCell(std::string);
    void VTKLineCell(std::string,std::string, std::string);
    void VTKLineCell(std::pair<std::string,std::string>,std::string);
    void VTKLineCell_set_background(std::string,std::string,std::string,_vec<double>);

    _vec<double> value_to_RGB(double,double,double,std::string);

    void VTKCell(Organ*);
    void VTKLine(Organ*);
    void VTKLineCell(Organ*);
    void VTK_example(void);
    void VTKLineCell_show_cell_number(Organ*,int,int);
}

#endif