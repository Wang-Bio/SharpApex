#include "../include/vtk2dv.h"

struct RGB{
    double r,g,b;
    
    RGB(double r=0, double g=0, double b=0) : r(r), g(g), b(b) {}
};

double lerp(double a, double b ,double t){
    return (1.0-t)*a + t*b;   
}

RGB palette_10_13_33(double t){
    RGB red(255,0,0);
    RGB light_red(250,100,0);
    RGB orange(255,185,0);
    RGB light_orange(210,235,65);
    RGB yellow(180,255,130);
    RGB light_green(140,235,190);
    RGB green(100,180,255);
    RGB light_blue(50,90,255);
    RGB blue(0,0,255);
    if(t<0){
        return red;
    }
    else if (t < 0.125) {
        double local_t = t / 0.125;
        return RGB(lerp(red.r, light_red.r, local_t), lerp(red.g, light_red.g, local_t), lerp(red.b, light_red.b, local_t));
    }
    else if(t < 0.25) {
        double local_t = (t - 0.125) / 0.125;
        return RGB(lerp(light_red.r, orange.r, local_t), lerp(light_red.g, orange.g, local_t), lerp(light_red.b, orange.b, local_t));
    }else if(t < 0.375){
        double local_t = (t - 0.25) / 0.125;
        return RGB(lerp(orange.r, light_orange.r, local_t), lerp(orange.g, light_orange.g, local_t), lerp(orange.b, light_orange.b, local_t));
    }
    else if (t < 0.5) {
        double local_t = (t - 0.375) / 0.125;
        return RGB(lerp(light_orange.r, yellow.r, local_t), lerp(light_orange.g, yellow.g, local_t), lerp(light_orange.b, yellow.b, local_t));
    }
    else if( t < 0.625){
        double local_t = (t - 0.5) / 0.125;
        return RGB(lerp(yellow.r, light_green.r, local_t), lerp(yellow.g, light_green.g, local_t), lerp(yellow.b, light_green.b, local_t));
    } else if (t < 0.75) {
        double local_t = (t - 0.625) / 0.125;
        return RGB(lerp(light_green.r, green.r, local_t), lerp(light_green.g, green.g, local_t), lerp(light_green.b, green.b, local_t));
    } 
    else if (t < 0.875) {
        double local_t = (t - 0.75) / 0.125;
        return RGB(lerp(green.r, light_blue.r, local_t), lerp(green.g, light_blue.g, local_t), lerp(green.b, light_blue.b, local_t));
    }
    else if (t < 1.0) {
        double local_t = (t - 0.875) / 0.125;
        return RGB(lerp(light_blue.r, blue.r, local_t), lerp(light_blue.g, blue.g, local_t), lerp(light_blue.b, blue.b, local_t));
    }
    else {
        return blue;
    }
}

namespace autoVTK{
    void VTKLine(std::string inputFile_line){
        //std::string inputFile_line = "/mnt/e/c_file/2023_06_20_biregion_angles_position_5000/y_0.60/y_0.60_ab_20.00/1/2dv_line0017900000.vtk";

        //Read the source file
        vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
        reader->SetFileName(inputFile_line.c_str());
        reader->Update();
        reader->GetOutput()->GetScalarRange();

        //create the mapper
        vtkSmartPointer<vtkDataSetMapper> mapper = vtkSmartPointer<vtkDataSetMapper>::New();
        mapper->SetInputConnection(reader->GetOutputPort());
        mapper->ScalarVisibilityOff(); //ignore scalar data

        //create the actor
        vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
        actor->SetMapper(mapper);
        actor->GetProperty()->SetColor(1.0,1.0,1.0); //set color to white
        actor->GetProperty()->SetLineWidth(0.1);

        //create the renderer
        vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
        renderer->AddActor(actor);
        renderer->SetBackground(0.32, 0.58, 0.67);

        //create the render window
        vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
        renderWindow->AddRenderer(renderer);

        //create the render window interactor
        vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
        renderWindowInteractor->SetRenderWindow(renderWindow);

        //render and interact
        renderWindow->Render();
        renderWindowInteractor->Start();
        //return EXIT_SUCCESS;
    }
    
    void VTKCell(std::string inputFile_cell){
        //std::string inputFile_cell = "/mnt/e/c_file/2023_06_20_biregion_angles_position_5000/y_0.60/y_0.60_ab_20.00/1/2dv_cell0017900000.vtk";

        //Read the source file
        vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
        reader->SetFileName(inputFile_cell.c_str());
        reader->Update();

        double range[2];
        reader->GetOutput()->GetScalarRange(range);

        //creat the color map
        vtkSmartPointer<vtkColorTransferFunction> colorFunc = vtkSmartPointer<vtkColorTransferFunction>::New();
        colorFunc->SetColorSpaceToDiverging();
        colorFunc->AddRGBPoint(range[0],0.230,0.299,0.754);
        colorFunc->AddRGBPoint(range[1],0.706,0.016,0.150);

        //create the mapper
        vtkSmartPointer<vtkDataSetMapper> mapper = vtkSmartPointer<vtkDataSetMapper>::New();
        mapper->SetInputConnection(reader->GetOutputPort());
        mapper->SetLookupTable(colorFunc);
        mapper->UseLookupTableScalarRangeOn();
        //mapper->ScalarVisibilityOff(); //ignore scalar data

        //create the actor
        vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
        actor->SetMapper(mapper);
        actor->GetProperty()->SetColor(1.0,1.0,1.0); //set color to white
        //actor->GetProperty()->SetLineWidth(0.1);

        //create the renderer
        vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
        renderer->AddActor(actor);
        renderer->SetBackground(0.32, 0.58, 0.67);

        //create the render window
        vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
        renderWindow->AddRenderer(renderer);

        //create the render window interactor
        vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
        renderWindowInteractor->SetRenderWindow(renderWindow);

        //render and interact
        renderWindow->Render();
        renderWindowInteractor->Start();
        //return EXIT_SUCCESS;
    }
    
    void VTKLineCell(std::string inputFile_line, std::string inputFile_cell, std::string output_png){
        //std::string inputFile_line = "/mnt/e/c_file/2023_06_20_biregion_angles_position_5000/y_0.60/y_0.60_ab_20.00/1/2dv_line0017900000.vtk";

        //Read the source file from line vtk
        vtkSmartPointer<vtkUnstructuredGridReader> reader_line = vtkSmartPointer<vtkUnstructuredGridReader>::New();
        reader_line->SetFileName(inputFile_line.c_str());
        reader_line->Update();
        reader_line->GetOutput()->GetScalarRange();

        //std::string inputFile_cell = "/mnt/e/c_file/2023_06_20_biregion_angles_position_5000/y_0.60/y_0.60_ab_20.00/1/2dv_cell0017900000.vtk";
        //Read the source file from cell vtk
        vtkSmartPointer<vtkUnstructuredGridReader> reader_cell = vtkSmartPointer<vtkUnstructuredGridReader>::New();
        reader_cell->SetFileName(inputFile_cell.c_str());
        reader_cell->Update();
        double range[2];
        reader_cell->GetOutput()->GetScalarRange(range);

        //creat the color map for cell vtk
        vtkSmartPointer<vtkColorTransferFunction> colorFunc = vtkSmartPointer<vtkColorTransferFunction>::New();
        colorFunc->SetColorSpaceToDiverging();
        colorFunc->AddRGBPoint(range[0],0.230,0.299,0.754);
        colorFunc->AddRGBPoint(range[1],0.706,0.016,0.150);

        //create the mapper for line and cell vtk
        vtkSmartPointer<vtkDataSetMapper> mapper_line = vtkSmartPointer<vtkDataSetMapper>::New();
        mapper_line->SetInputConnection(reader_line->GetOutputPort());
        mapper_line->ScalarVisibilityOff(); //ignore scalar data

        vtkSmartPointer<vtkDataSetMapper> mapper_cell = vtkSmartPointer<vtkDataSetMapper>::New();
        mapper_cell->SetInputConnection(reader_cell->GetOutputPort());
        mapper_cell->SetLookupTable(colorFunc);
        mapper_cell->UseLookupTableScalarRangeOn();
        //mapper->ScalarVisibilityOff(); //ignore scalar data
        
        //create the actor
        vtkSmartPointer<vtkActor> actor_line = vtkSmartPointer<vtkActor>::New();
        actor_line->SetMapper(mapper_line);
        actor_line->GetProperty()->SetColor(1.0,1.0,1.0); //set color to white
        actor_line->GetProperty()->SetLineWidth(0.1);

        vtkSmartPointer<vtkActor> actor_cell = vtkSmartPointer<vtkActor>::New();
        actor_cell->SetMapper(mapper_cell);
        actor_cell->GetProperty()->SetColor(1.0,1.0,1.0); //set color to white
        //actor->GetProperty()->SetLineWidth(0.1);

        //create the renderer
        vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
        renderer->AddActor(actor_cell);
        renderer->AddActor(actor_line);
        // Enable a gradient background
        renderer->GradientBackgroundOn();

        // Set the background colors
        // The top color (light blue) is RGB(0.32, 0.58, 0.67)
        renderer->SetBackground(0.32, 0.58, 0.67);
        // The bottom color (dark blue) is RGB(0.20, 0.41, 0.51)
        renderer->SetBackground2(0.20, 0.41, 0.51);


        //create a light
        vtkSmartPointer<vtkLight> light = vtkSmartPointer<vtkLight>::New();
        light->SetIntensity(0.9);
        //light->SetAmbientColor(0.5,0.5,0.5);
        renderer->AddLight(light);

        //create the render window
        vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
        renderWindow->AddRenderer(renderer);
        renderWindow->SetSize(1080,1080);

        //create the render window interactor
        vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
        renderWindowInteractor->SetRenderWindow(renderWindow);

        //render and interact
        renderWindowInteractor->Initialize();
        renderWindow->OffScreenRenderingOn();
        renderWindow->Render();
        //renderWindowInteractor->Start();


        // Capture the image
        vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter =
            vtkSmartPointer<vtkWindowToImageFilter>::New();
        windowToImageFilter->SetInput(renderWindow);
        windowToImageFilter->ReadFrontBufferOff();
        windowToImageFilter->Update();

        // Write the image to a file
        vtkSmartPointer<vtkPNGWriter> writer =
            vtkSmartPointer<vtkPNGWriter>::New();
        writer->SetFileName(output_png.c_str());
        writer->SetInputConnection(windowToImageFilter->GetOutputPort());
        writer->Write();
        //return EXIT_SUCCESS;
    }

    void VTKLineCell(std::pair<std::string,std::string> inputFile_cell_line, std::string output_png){
        //std::string inputFile_line = "/mnt/e/c_file/2023_06_20_biregion_angles_position_5000/y_0.60/y_0.60_ab_20.00/1/2dv_line0017900000.vtk";

        //Read the source file from line vtk
        vtkSmartPointer<vtkUnstructuredGridReader> reader_line = vtkSmartPointer<vtkUnstructuredGridReader>::New();
        reader_line->SetFileName(inputFile_cell_line.second.c_str());
        reader_line->Update();
        reader_line->GetOutput()->GetScalarRange();

        //std::string inputFile_cell = "/mnt/e/c_file/2023_06_20_biregion_angles_position_5000/y_0.60/y_0.60_ab_20.00/1/2dv_cell0017900000.vtk";
        //Read the source file from cell vtk
        vtkSmartPointer<vtkUnstructuredGridReader> reader_cell = vtkSmartPointer<vtkUnstructuredGridReader>::New();
        reader_cell->SetFileName(inputFile_cell_line.first.c_str());
        reader_cell->Update();
        double range[2];
        reader_cell->GetOutput()->GetScalarRange(range);

        //creat the color map for cell vtk
        vtkSmartPointer<vtkColorTransferFunction> colorFunc = vtkSmartPointer<vtkColorTransferFunction>::New();
        colorFunc->SetColorSpaceToDiverging();
        colorFunc->AddRGBPoint(range[0],0.230,0.299,0.754);
        colorFunc->AddRGBPoint(range[1],0.706,0.016,0.150);

        //create the mapper for line and cell vtk
        vtkSmartPointer<vtkDataSetMapper> mapper_line = vtkSmartPointer<vtkDataSetMapper>::New();
        mapper_line->SetInputConnection(reader_line->GetOutputPort());
        mapper_line->ScalarVisibilityOff(); //ignore scalar data

        vtkSmartPointer<vtkDataSetMapper> mapper_cell = vtkSmartPointer<vtkDataSetMapper>::New();
        mapper_cell->SetInputConnection(reader_cell->GetOutputPort());
        mapper_cell->SetLookupTable(colorFunc);
        mapper_cell->UseLookupTableScalarRangeOn();
        //mapper->ScalarVisibilityOff(); //ignore scalar data
        
        //create the actor
        vtkSmartPointer<vtkActor> actor_line = vtkSmartPointer<vtkActor>::New();
        actor_line->SetMapper(mapper_line);
        actor_line->GetProperty()->SetColor(1.0,1.0,1.0); //set color to white
        actor_line->GetProperty()->SetLineWidth(0.1);

        vtkSmartPointer<vtkActor> actor_cell = vtkSmartPointer<vtkActor>::New();
        actor_cell->SetMapper(mapper_cell);
        actor_cell->GetProperty()->SetColor(1.0,1.0,1.0); //set color to white
        //actor->GetProperty()->SetLineWidth(0.1);

        //create the renderer
        vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
        renderer->AddActor(actor_cell);
        renderer->AddActor(actor_line);
        // Enable a gradient background
        renderer->GradientBackgroundOn();

        // Set the background colors
        // The top color (light blue) is RGB(0.32, 0.58, 0.67)
        renderer->SetBackground(0.32, 0.58, 0.67);
        // The bottom color (dark blue) is RGB(0.20, 0.41, 0.51)
        renderer->SetBackground2(0.20, 0.41, 0.51);


        //create a light
        vtkSmartPointer<vtkLight> light = vtkSmartPointer<vtkLight>::New();
        light->SetIntensity(0.9);
        //light->SetAmbientColor(0.5,0.5,0.5);
        renderer->AddLight(light);

        //create the render window
        vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
        renderWindow->AddRenderer(renderer);
        renderWindow->SetSize(1080,1080);

        //create the render window interactor
        vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
        renderWindowInteractor->SetRenderWindow(renderWindow);

        //render and interact
        renderWindowInteractor->Initialize();
        renderWindow->OffScreenRenderingOn();
        renderWindow->Render();
        //renderWindowInteractor->Start();


        // Capture the image
        vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter =
            vtkSmartPointer<vtkWindowToImageFilter>::New();
        windowToImageFilter->SetInput(renderWindow);
        windowToImageFilter->ReadFrontBufferOff();
        windowToImageFilter->Update();

        // Write the image to a file
        vtkSmartPointer<vtkPNGWriter> writer =
            vtkSmartPointer<vtkPNGWriter>::New();
        writer->SetFileName(output_png.c_str());
        writer->SetInputConnection(windowToImageFilter->GetOutputPort());
        writer->Write();
        //return EXIT_SUCCESS;
    }

    void VTKLineCell_set_background(std::string inputFile_line, std::string inputFile_cell, std::string output_png, _vec<double> RGB_background){
        //std::string inputFile_line = "/mnt/e/c_file/2023_06_20_biregion_angles_position_5000/y_0.60/y_0.60_ab_20.00/1/2dv_line0017900000.vtk";

        //Read the source file from line vtk
        vtkSmartPointer<vtkUnstructuredGridReader> reader_line = vtkSmartPointer<vtkUnstructuredGridReader>::New();
        reader_line->SetFileName(inputFile_line.c_str());
        reader_line->Update();
        reader_line->GetOutput()->GetScalarRange();

        //std::string inputFile_cell = "/mnt/e/c_file/2023_06_20_biregion_angles_position_5000/y_0.60/y_0.60_ab_20.00/1/2dv_cell0017900000.vtk";
        //Read the source file from cell vtk
        vtkSmartPointer<vtkUnstructuredGridReader> reader_cell = vtkSmartPointer<vtkUnstructuredGridReader>::New();
        reader_cell->SetFileName(inputFile_cell.c_str());
        reader_cell->Update();
        double range[2];
        reader_cell->GetOutput()->GetScalarRange(range);

        //creat the color map for cell vtk
        vtkSmartPointer<vtkColorTransferFunction> colorFunc = vtkSmartPointer<vtkColorTransferFunction>::New();
        colorFunc->SetColorSpaceToDiverging();
        colorFunc->AddRGBPoint(range[0],0.230,0.299,0.754);
        colorFunc->AddRGBPoint(range[1],0.706,0.016,0.150);

        //create the mapper for line and cell vtk
        vtkSmartPointer<vtkDataSetMapper> mapper_line = vtkSmartPointer<vtkDataSetMapper>::New();
        mapper_line->SetInputConnection(reader_line->GetOutputPort());
        mapper_line->ScalarVisibilityOff(); //ignore scalar data

        vtkSmartPointer<vtkDataSetMapper> mapper_cell = vtkSmartPointer<vtkDataSetMapper>::New();
        mapper_cell->SetInputConnection(reader_cell->GetOutputPort());
        mapper_cell->SetLookupTable(colorFunc);
        mapper_cell->UseLookupTableScalarRangeOn();
        //mapper->ScalarVisibilityOff(); //ignore scalar data
        
        //create the actor
        vtkSmartPointer<vtkActor> actor_line = vtkSmartPointer<vtkActor>::New();
        actor_line->SetMapper(mapper_line);
        actor_line->GetProperty()->SetColor(1.0,1.0,1.0); //set color to white
        actor_line->GetProperty()->SetLineWidth(0.1);

        vtkSmartPointer<vtkActor> actor_cell = vtkSmartPointer<vtkActor>::New();
        actor_cell->SetMapper(mapper_cell);
        actor_cell->GetProperty()->SetColor(1.0,1.0,1.0); //set color to white
        //actor->GetProperty()->SetLineWidth(0.1);

        //create the renderer
        vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
        renderer->AddActor(actor_cell);
        renderer->AddActor(actor_line);
        // Enable a gradient background
        //renderer->GradientBackgroundOn();

        // Set the background colors
        // The top color (light blue) is RGB(0.32, 0.58, 0.67)
        renderer->SetBackground(RGB_background.x/255, RGB_background.y/255, RGB_background.z/255);



        //create a light
        vtkSmartPointer<vtkLight> light = vtkSmartPointer<vtkLight>::New();
        light->SetIntensity(0.9);
        light->SetAmbientColor(0.5,0.5,0.5);
        renderer->AddLight(light);

        //create the render window
        vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
        renderWindow->AddRenderer(renderer);
        renderWindow->SetSize(1080,1080);

        //create the render window interactor
        vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
        renderWindowInteractor->SetRenderWindow(renderWindow);

        //render and interact
        renderWindowInteractor->Initialize();
        renderWindow->OffScreenRenderingOn();
        renderWindow->Render();
        //renderWindowInteractor->Start();


        // Capture the image
        vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter =
            vtkSmartPointer<vtkWindowToImageFilter>::New();
        windowToImageFilter->SetInput(renderWindow);
        windowToImageFilter->ReadFrontBufferOff();
        windowToImageFilter->Update();

        // Write the image to a file
        vtkSmartPointer<vtkPNGWriter> writer =
            vtkSmartPointer<vtkPNGWriter>::New();
        writer->SetFileName(output_png.c_str());
        writer->SetInputConnection(windowToImageFilter->GetOutputPort());
        writer->Write();
        //return EXIT_SUCCESS;
    }

    _vec<double> value_to_RGB(double value, double maximum_value, double minimum_value, std::string palette_mode){
        _vec<double> RGB_value;
        RGB RGB_value_;
        double relative_value = (value-minimum_value)/(maximum_value-minimum_value);
        if(palette_mode == "10,13,33"){
            RGB_value_ = palette_10_13_33(relative_value);
        }
        RGB_value.x = RGB_value_.r;
        RGB_value.y = RGB_value_.g;
        RGB_value.z = RGB_value_.b;
        cout<<"relative_value: "<<relative_value<<"r: "<< RGB_value.x<<"; g: "<<RGB_value.y<<"; b: "<<RGB_value.z<<std::endl;
        return RGB_value;
    }

    void VTK_example(void){
            vtkSmartPointer<vtkConeSource> coneSource = vtkSmartPointer<vtkConeSource>::New();
            coneSource->Update();

            // Map the cone to an actor
            vtkSmartPointer<vtkPolyDataMapper> coneMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
            coneMapper->SetInputConnection(coneSource->GetOutputPort());

            vtkSmartPointer<vtkActor> coneActor = vtkSmartPointer<vtkActor>::New();
            coneActor->SetMapper(coneMapper);

            // Create a renderer
            vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
            renderer->AddActor(coneActor);

            // Create a render window and render window interactor
            vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
            renderWindow->AddRenderer(renderer);

            vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
            renderWindowInteractor->SetRenderWindow(renderWindow);

            // Start the interaction
            renderWindow->Render();
            renderWindowInteractor->Start();

    }

    void VTKCell(Organ* p_g){
        // Create points object
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

        // Create cell array to store polygons
        vtkSmartPointer<vtkCellArray> polygons = vtkSmartPointer<vtkCellArray>::New();

        // The points and polygons are added here
        for(int vi=0; vi<(int)p_g->p_v.size(); vi++){
            points->InsertNextPoint(p_g->p_v[vi]->loc.x, p_g->p_v[vi]->loc.y, p_g->p_v[vi]->loc.z);
        }
        for(int ci=0; ci<(int)p_g->p_c.size(); ci++){
            polygons->InsertNextCell(p_g->p_c[ci]->vi.size());
            for(int vi=0; vi<(int)p_g->p_c[ci]->vi.size(); vi++){
                polygons->InsertCellPoint(p_g->p_c[ci]->vi[vi]);
            }
        }

        // Create vtkDoubleArray for scalar values (served for color of polygons)
        std::vector<int> scalarValues;

        vtkSmartPointer<vtkIntArray> scalars = vtkSmartPointer<vtkIntArray>::New();
        for(int ci=0; ci<(int)p_g->p_c.size(); ci++){
            scalarValues.push_back(p_g->p_c[ci]->cellDivisionCount);
            scalars->InsertNextValue(p_g->p_c[ci]->cellDivisionCount);
        }

        // Create polydata object and set points and polygons
        vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
        polyData->SetPoints(points);
        polyData->SetPolys(polygons);
        //Assign the scalars to the polyData
        polyData->GetCellData()->SetScalars(scalars);

        //Create a lookup table to map scalar values to colors
        int minValue = *std::min_element(scalarValues.begin(), scalarValues.end());
        int maxValue = *std::max_element(scalarValues.begin(), scalarValues.end());

        vtkSmartPointer<vtkLookupTable> lut = vtkSmartPointer<vtkLookupTable>::New();
        lut->SetRange(minValue, maxValue); 
        lut->Build(); 

        // Create mapper
        vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper->SetInputData(polyData);
        
        // Link the lookup table with the mapper
        mapper->SetLookupTable(lut);
        mapper->SetScalarRange(minValue, maxValue);
        mapper->SetScalarModeToUseCellData();


        // Create actor
        vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
        actor->SetMapper(mapper);

        // Create renderer
        vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();

        // Create render window
        vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
        renderWindow->AddRenderer(renderer);

        // Create render window interactor
        vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
        renderWindowInteractor->SetRenderWindow(renderWindow);

        // Add actor to renderer
        renderer->AddActor(actor);
        renderer->SetBackground(0.32, 0.58, 0.67);

        // Start interaction
        renderWindow->Render();
        renderWindowInteractor->Start();
    }

    void VTKLine(Organ* p_g){
        // Create points object
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

        // Create cell array to store lines
        vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();

        // The points and lines are added here
        for(int vi=0; vi<(int)p_g->p_v.size(); vi++){
            points->InsertNextPoint(p_g->p_v[vi]->loc.x, p_g->p_v[vi]->loc.y, p_g->p_v[vi]->loc.z);
        }
        for(int li=0; li<(int)p_g->p_l.size(); li++){
            lines->InsertNextCell(2);
            //std::cout<<"p_g->p_l[li]->vi[0]: "<<p_g->p_l[li]->vi[0]<<",p_g->p_l[li]->vi[1]: "<<p_g->p_l[li]->vi[1]<<std::endl;
            lines->InsertCellPoint(p_g->p_l[li]->vi[0]);
            lines->InsertCellPoint(p_g->p_l[li]->vi[1]);
        }

        // Create polydata object and set points and polygons
        vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
        polyData->SetPoints(points);
        polyData->SetLines(lines);

        // Create mapper
        vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper->SetInputData(polyData);

        // Create actor
        vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
        actor->SetMapper(mapper);
        actor->GetProperty()->SetColor(1.0,1.0,1.0); //set color to white
        actor->GetProperty()->SetLineWidth(1.0);
        // Create renderer
        vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();

        // Create render window
        vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
        renderWindow->AddRenderer(renderer);

        // Create render window interactor
        vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
        renderWindowInteractor->SetRenderWindow(renderWindow);

        // Add actor to renderer
        renderer->AddActor(actor);
        renderer->SetBackground(0.32, 0.58, 0.67);

        // Start interaction
        renderWindow->Render();
        renderWindowInteractor->Start();
    }

    void VTKLineCell(Organ* p_g){
    //1. collect coordinates of all points, and connectivity of all polygons and lines
        // Create points object
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

        // Create cell array to store polygons
        vtkSmartPointer<vtkCellArray> polygons = vtkSmartPointer<vtkCellArray>::New();

        // Create cell array to store lines
        vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
        
        // The points are added here
        for(int vi=0; vi<(int)p_g->p_v.size(); vi++){
            points->InsertNextPoint(p_g->p_v[vi]->loc.x, p_g->p_v[vi]->loc.y, p_g->p_v[vi]->loc.z);
        }
        // The polygons are added here
        for(int ci=0; ci<(int)p_g->p_c.size(); ci++){
            polygons->InsertNextCell(p_g->p_c[ci]->vi.size());
            for(int vi=0; vi<(int)p_g->p_c[ci]->vi.size(); vi++){
                polygons->InsertCellPoint(p_g->p_c[ci]->vi[vi]);
            }
        }
        // The lines are added here
        for(int li=0; li<(int)p_g->p_l.size(); li++){
            lines->InsertNextCell(2);
            lines->InsertCellPoint(p_g->p_l[li]->vi[0]);
            lines->InsertCellPoint(p_g->p_l[li]->vi[1]);
        }
        // Create polydata object and set points, lines, and polygons
        vtkSmartPointer<vtkPolyData> polyData_cell = vtkSmartPointer<vtkPolyData>::New();
        vtkSmartPointer<vtkPolyData> polyData_line = vtkSmartPointer<vtkPolyData>::New();

        polyData_cell->SetPoints(points);
        polyData_cell->SetPolys(polygons);
        polyData_line->SetPoints(points);
        polyData_line->SetLines(lines);
    //2. color scale for cells (lookup table)
        // Create vtkDoubleArray for scalar values (served for color of polygons)
        std::vector<int> scalarValues;
        vtkSmartPointer<vtkIntArray> scalars = vtkSmartPointer<vtkIntArray>::New();
        for(int ci=0; ci<(int)p_g->p_c.size(); ci++){
            scalarValues.push_back(p_g->p_c[ci]->cellDivisionCount);
            scalars->InsertNextValue(p_g->p_c[ci]->cellDivisionCount);
        }
        //Create a lookup table to map scalar values to colors
        int minValue = *std::min_element(scalarValues.begin(), scalarValues.end());
        int maxValue = *std::max_element(scalarValues.begin(), scalarValues.end());

        vtkSmartPointer<vtkLookupTable> lut = vtkSmartPointer<vtkLookupTable>::New();
        lut->SetRange(minValue, maxValue); 
        lut->Build();

        //create mapper
        vtkSmartPointer<vtkPolyDataMapper> mapper_line = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper_line->SetInputData(polyData_line);
        vtkSmartPointer<vtkPolyDataMapper> mapper_cell = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper_cell->SetInputData(polyData_cell);
        
        // Link the lookup table with the mapper
        mapper_cell->SetLookupTable(lut);
        //mapper_cell->SetScalarRange(minValue, maxValue);
        mapper_cell->SetScalarRange(0, 10);
        mapper_cell->SetScalarModeToUseCellData();

        // Create actor
        vtkSmartPointer<vtkActor> actor_cell = vtkSmartPointer<vtkActor>::New();
        vtkSmartPointer<vtkActor> actor_line = vtkSmartPointer<vtkActor>::New();
        actor_cell->SetMapper(mapper_cell);
        actor_line->SetMapper(mapper_line);
        actor_line->GetProperty()->SetColor(1.0,1.0,1.0); //set color to white
        actor_line->GetProperty()->SetLineWidth(1.0);

        // Create renderer
        vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();

        // Create render window
        vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
        renderWindow->AddRenderer(renderer);

        // Create render window interactor
        vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
        renderWindowInteractor->SetRenderWindow(renderWindow);

        // Add actor to renderer
        //renderer->AddActor(actor_line);
        renderer->AddActor(actor_cell);
        renderer->SetBackground(0.32, 0.58, 0.67);

        // Start interaction
        renderWindow->Render();
        renderWindowInteractor->Start();
        // Wait for 2 seconds
        //vtkTimerLog::WaitForSeconds(2);

        // Close the render window
        //renderWindowInteractor->GetRenderWindow()->Finalize();
    }

    void VTKLineCell_show_cell_number(Organ* p_g, int min_cbvalue, int max_cbvalue){
        //1. collect coordinates of all points, and connectivity of all polygons and lines
            // Create points object
            vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

            // Create cell array to store polygons
            vtkSmartPointer<vtkCellArray> polygons = vtkSmartPointer<vtkCellArray>::New();

            // Create cell array to store lines
            vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
            
            // The points are added here
            for(int vi=0; vi<(int)p_g->p_v.size(); vi++){
                points->InsertNextPoint(p_g->p_v[vi]->loc.x, p_g->p_v[vi]->loc.y, p_g->p_v[vi]->loc.z);
            }
            // The polygons are added here
            for(int ci=0; ci<(int)p_g->p_c.size(); ci++){
                polygons->InsertNextCell(p_g->p_c[ci]->vi.size());
                for(int vi=0; vi<(int)p_g->p_c[ci]->vi.size(); vi++){
                    polygons->InsertCellPoint(p_g->p_c[ci]->vi[vi]);
                }
            }
            // The lines are added here
            for(int li=0; li<(int)p_g->p_l.size(); li++){
                lines->InsertNextCell(2);
                lines->InsertCellPoint(p_g->p_l[li]->vi[0]);
                lines->InsertCellPoint(p_g->p_l[li]->vi[1]);
            }
            // Create polydata object and set points, lines, and polygons
            vtkSmartPointer<vtkPolyData> polyData_cell = vtkSmartPointer<vtkPolyData>::New();
            vtkSmartPointer<vtkPolyData> polyData_line = vtkSmartPointer<vtkPolyData>::New();

            polyData_cell->SetPoints(points);
            polyData_cell->SetPolys(polygons);
            polyData_line->SetPoints(points);
            polyData_line->SetLines(lines);
        //2. color scale for cells (lookup table)
            // Create vtkDoubleArray for scalar values (served for color of polygons)
            std::vector<int> scalarValues;
            vtkSmartPointer<vtkIntArray> scalars = vtkSmartPointer<vtkIntArray>::New();
            for(int ci=0; ci<(int)p_g->p_c.size(); ci++){
                scalarValues.push_back(p_g->p_c[ci]->cellDivisionCount);
                scalars->InsertNextValue(p_g->p_c[ci]->cellDivisionCount);
            }
            //Create a lookup table to map scalar values to colors
            int minValue = min_cbvalue;
            int maxValue = max_cbvalue;

            vtkSmartPointer<vtkLookupTable> lut = vtkSmartPointer<vtkLookupTable>::New();
            lut->SetRange(minValue, maxValue); 
            lut->Build();

            //create mapper
            vtkSmartPointer<vtkPolyDataMapper> mapper_line = vtkSmartPointer<vtkPolyDataMapper>::New();
            mapper_line->SetInputData(polyData_line);
            vtkSmartPointer<vtkPolyDataMapper> mapper_cell = vtkSmartPointer<vtkPolyDataMapper>::New();
            mapper_cell->SetInputData(polyData_cell);
            
            // Link the lookup table with the mapper
            mapper_cell->SetLookupTable(lut);
            //mapper_cell->SetScalarRange(minValue, maxValue);
            mapper_cell->SetScalarRange(0, 10);
            mapper_cell->SetScalarModeToUseCellData();

            // Create actor
            vtkSmartPointer<vtkActor> actor_cell = vtkSmartPointer<vtkActor>::New();
            vtkSmartPointer<vtkActor> actor_line = vtkSmartPointer<vtkActor>::New();
            actor_cell->SetMapper(mapper_cell);
            actor_line->SetMapper(mapper_line);
            actor_line->GetProperty()->SetColor(1.0,1.0,1.0); //set color to white
            actor_line->GetProperty()->SetLineWidth(1.0);

            // Create renderer
            vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();

            // Create render window
            vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
            renderWindow->AddRenderer(renderer);

            // Create render window interactor
            vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
            renderWindowInteractor->SetRenderWindow(renderWindow);

            // Add actor to renderer
            //renderer->AddActor(actor_line);
            renderer->AddActor(actor_cell);
            renderer->SetBackground(0.32, 0.58, 0.67);

            // Create a text actor
            vtkSmartPointer<vtkTextActor> textActor = vtkSmartPointer<vtkTextActor>::New();
            char cell_number_vtk[100];
            sprintf(cell_number_vtk,"Cell number : %d",(int)p_g->p_c.size());
            textActor->SetInput(cell_number_vtk);
            textActor->GetTextProperty()->SetFontSize(24); // Adjust as needed
            textActor->GetTextProperty()->SetColor(1.0, 1.0, 1.0); // Color: White

            // Set the position of the text actor. (0, 0) is bottom-left, (1, 1) is top-right
            textActor->SetPosition(0.5, 0.95); // Near top-middle corner

            // Add the text actor to the renderer
            renderer->AddActor(textActor);
            // Start interaction
            renderWindow->Render();
            renderWindowInteractor->Start();
            // Wait for 2 seconds


            // Close the render window
            renderWindow->Finalize();
        }

}