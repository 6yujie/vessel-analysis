// VTK: Spiral with vtkTubeFilter
// Varying tube radius and independent RGB colors with an unsignedCharArray
// Contributed by Marcus Thamson

#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>

#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkDataSet.h>
#include <vtkDataSetAttributes.h>
#include <vtkProperty.h>
#include <vtkSmartPointer.h>
#include <vtkTubeFilter.h>

#include <vtkTriangleStrip.h>
#include <vtkTriangle.h>

#include <vtkPolyDataNormals.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>

#include <vtkDataSetMapper.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkCamera.h>
#include <vtkInteractorStyleTrackballCamera.h>

#include <vtkWindowToImageFilter.h>
#include <vtkUnsignedCharArray.h>

#include <QVTKOpenGLWidget.h>

#include <vtkMath.h>

#include <QApplication>
#include <QMainWindow>
#include <QSurfaceFormat>
#include <QVector3D>
#include <QVector>
#include <QDebug>
#include <QImage>

#include "glwidget.h"

int main(int argc, char **argv)
{
	// Spiral tube
	double vX, vY, vZ;
// 	unsigned int nV = 256;      // No. of vertices
	unsigned int nV = 50;      // No. of vertices
// 	unsigned int nCyc = 5;      // No. of spiral cycles
	unsigned int nCyc = 1;      // No. of spiral cycles
	double rT1 = 0.1, rT2 = 0.5;// Start/end tube radii
	double rS = 2;              // Spiral radius
	double h = 10;              // Height
	unsigned int nTv = 8;       // No. of surface elements for each tube vertex

	unsigned int i;

	// Create points and cells for the spiral
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	for (i = 0; i < nV; i++)
	{
		// Spiral coordinates
		vX = rS * cos(2 * vtkMath::Pi() * nCyc * i / (nV - 1));
		vY = rS * sin(2 * vtkMath::Pi() * nCyc * i / (nV - 1));
		vZ = h * i / nV;
		points->InsertPoint(i, vX, vY, vZ);
	}

	// Cell: 单元，一系列有序的点按照指定类型连接所定义的结构
	vtkSmartPointer<vtkCellArray> lines =
		vtkSmartPointer<vtkCellArray>::New();
	lines->InsertNextCell(nV);
	for (i = 0; i < nV; i++)
	{
		lines->InsertCellPoint(i);
	}

	// vtkPolyData 表示由顶点、线、多边形和/或三角形组成的几何结构
	vtkSmartPointer<vtkPolyData> polyData =
		vtkSmartPointer<vtkPolyData>::New();
	polyData->SetPoints(points);
	polyData->SetLines(lines);

	// Varying tube radius using sine-function
	vtkSmartPointer<vtkDoubleArray> tubeRadius =
		vtkSmartPointer<vtkDoubleArray>::New();
	tubeRadius->SetName("TubeRadius");
	tubeRadius->SetNumberOfTuples(nV);
	for (i = 0; i < nV; i++)
	{
		tubeRadius->SetTuple1(i,
			rT1 + (rT2 - rT1) * sin(vtkMath::Pi() * i / (nV - 1)));
	}
	polyData->GetPointData()->AddArray(tubeRadius);
	polyData->GetPointData()->SetActiveScalars("TubeRadius");

	// RBG array (could add Alpha channel too I guess...)
	// Varying from blue to red
	vtkSmartPointer<vtkUnsignedCharArray> colors =
		vtkSmartPointer<vtkUnsignedCharArray>::New();
	colors->SetName("Colors");
	colors->SetNumberOfComponents(3);
	colors->SetNumberOfTuples(nV);
	for (i = 0; i < nV; i++)
	{
		colors->InsertTuple3(i,
			int(255 * i / (nV - 1)),
			0,
			int(255 * (nV - 1 - i) / (nV - 1)));
	}
	polyData->GetPointData()->AddArray(colors);

	vtkSmartPointer<vtkTubeFilter> tube
		= vtkSmartPointer<vtkTubeFilter>::New();
	tube->SetInputData(polyData);
	tube->SetNumberOfSides(nTv);
	tube->SetVaryRadiusToVaryRadiusByAbsoluteScalar();
	tube->Update();

	vtkSmartPointer<vtkPolyDataMapper> mapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputConnection(tube->GetOutputPort());
	mapper->ScalarVisibilityOn();
	mapper->SetScalarModeToUsePointFieldData();
	mapper->SelectColorArray("Colors");

// 	vtkPolyDataNormals *normals = vtkPolyDataNormals::New();
	vtkPolyData* tbData = tube->GetOutput();

	vtkPoints *tbPoints = tbData->GetPoints();
	vtkCellArray *tbPolys = tbData->GetPolys();
	vtkCellArray *tbStrips = tbData->GetStrips();

	qDebug() << "points: " << tbPoints->GetNumberOfPoints(); // 400
	qDebug() << "polys: " << tbPolys->GetNumberOfCells(); // 0
	qDebug() << "strips: " << tbStrips->GetNumberOfCells(); // 8

	vtkIdType npts = 0;
	vtkIdType *indx = nullptr;
	vtkCellArray *tbPolyStrips = vtkCellArray::New();
	if (tbStrips->GetNumberOfCells() > 0)
	{
		vtkIdType *ptIds = nullptr;
		for (tbStrips->InitTraversal(); tbStrips->GetNextCell(npts, ptIds);)
		{
			vtkTriangleStrip::DecomposeStrip(npts, ptIds, tbPolyStrips);
		}
	}

	QVector<QVector3D> m_vertices;
	QVector<QVector3D> m_normals;

	double n[3], v1[3], v2[3], v3[3];
	for (tbPolyStrips->InitTraversal(); tbPolyStrips->GetNextCell(npts, indx);)
	{
		tbPoints->GetPoint(indx[0], v1);
		tbPoints->GetPoint(indx[1], v2);
		tbPoints->GetPoint(indx[2], v3);

		vtkTriangle::ComputeNormal(tbPoints, npts, indx, n);
		m_vertices << QVector3D(v1[0], v1[1], v1[2]);
		m_vertices << QVector3D(v2[0], v2[1], v2[2]);
		m_vertices << QVector3D(v3[0], v3[1], v3[2]);
		m_normals << QVector3D(n[0], n[1], n[2]);
		m_normals << QVector3D(n[0], n[1], n[2]);
		m_normals << QVector3D(n[0], n[1], n[2]);
	}

#if 0
	vtkDataArray* vertexes = tbPoints->GetData();
	vtkDataArray* normals = tbData->GetPointData()->GetNormals();

	auto vertexLen = vertexes->GetNumberOfValues();
	auto normalLen = normals->GetNumberOfValues();
	QVector<QVector3D> m_vertices;
	QVector<QVector3D> m_normals;

	if (vertexLen == normalLen)
	{
		for (int i = 0; i < vertexLen; i++)
		{
			float vx = vertexes->GetComponent(i, 0);
			float vy = vertexes->GetComponent(i, 1);
			float vz = vertexes->GetComponent(i, 2);
			m_vertices << QVector3D(vx, vy, vz);

			float nx = normals->GetComponent(i, 0);
			float ny = normals->GetComponent(i, 1);
			float nz = normals->GetComponent(i, 2);
			m_normals << QVector3D(nx, ny, nz);
// 			m_normals << QVector3D(0, 0, 0);
		}

// 		qDebug() << m_vertices;
		
	}
#endif

// 	normals->SetInputData(tube-());

	vtkSmartPointer<vtkActor> actor =
		vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);

	vtkSmartPointer<vtkRenderer> renderer =
		vtkSmartPointer<vtkRenderer>::New();
	renderer->AddActor(actor);
	renderer->SetBackground(.2, .3, .4);

	// Make an oblique view
	renderer->GetActiveCamera()->Azimuth(30);
	renderer->GetActiveCamera()->Elevation(30);
	renderer->ResetCamera();

	vtkSmartPointer<vtkRenderWindow> renWin =
		vtkSmartPointer<vtkRenderWindow>::New();
	vtkSmartPointer<vtkRenderWindowInteractor>
		iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();

	iren->SetRenderWindow(renWin);
	renWin->AddRenderer(renderer);
	renWin->SetSize(500, 500);
// 	renWin->OffScreenRenderingOn();
	renWin->Render();

	vtkSmartPointer<vtkInteractorStyleTrackballCamera> style =
		vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
	iren->SetInteractorStyle(style);

	vtkWindowToImageFilter *w2i = vtkWindowToImageFilter::New();
	renWin->SwapBuffersOff();
	w2i->SetInput(renWin);
// 	w2i->SetScale(2);
	w2i->ReadFrontBufferOff();
	w2i->Update();

// 	vtkPNGWriter *imgWriter = vtkPNGWriter::New();
// 	imgWriter->SetInputConnection(w2i->GetOutputPort());
// 	imgWriter->WriteToMemoryOn();

	iren->Start();
// 
// 	imgWriter->Update();
// 	imgWriter->Write();
// 	vtkUnsignedCharArray* result = imgWriter->GetResult();

#if 0
	vtkImageData* img = w2i->GetOutput();
	int dims[3];
	img->GetDimensions(dims);
	vtkUnsignedCharArray* result = vtkUnsignedCharArray::SafeDownCast(img->GetPointData()->GetScalars());
// 	vtkDataArray* vtk_array = img->GetPointData()->GetScalars();

	QImage qImg(dims[0], dims[1], QImage::Format_ARGB32);
	vtkIdType tupleIndex = 0;
	int qImageBitIndex = 0;
	QRgb *qImageBits = (QRgb*)qImg.bits();
	unsigned char* scalarTuples = result->GetPointer(0);

	for (int j = 0; j < dims[1]; j++)
	{
		for (int i = 0; i < dims[0]; i++)
		{
			unsigned char *tuple = scalarTuples + (tupleIndex++ * 3);
			QRgb color = qRgba(tuple[0], tuple[1], tuple[2], 255);
			*(qImageBits + (qImageBitIndex++)) = color;
		}
	}

// 	qImg.save("myImage.jpg");
#endif


	QApplication a(argc, argv);

	QMainWindow w;


	GLWidget *gw = new GLWidget(m_vertices, m_normals);

	w.setCentralWidget(gw);
	w.resize(600, 500);
	w.show();

	return a.exec();
}
