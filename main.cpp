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

#include <vtkPolyDataNormals.h>

#include <vtkDataSetMapper.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkCamera.h>
#include <vtkInteractorStyleTrackballCamera.h>

#include <vtkMath.h>

#include <QApplication>
#include <QMainWindow>
#include <QSurfaceFormat>
#include <QVector3D>
#include <QVector>
#include <QDebug>

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
// 			m_normals << QVector3D(nx, ny, nz);
			m_normals << QVector3D(0, 0, 0);
		}

// 		qDebug() << m_vertices;
		
	}

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
	renWin->Render();

	vtkSmartPointer<vtkInteractorStyleTrackballCamera> style =
		vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
	iren->SetInteractorStyle(style);

	iren->Start();

	QApplication a(argc, argv);

	QMainWindow w;

	GLWidget *gw = new GLWidget(m_vertices, m_normals);

	w.setCentralWidget(gw);
	w.resize(600, 500);
	w.show();

	return a.exec();
}
