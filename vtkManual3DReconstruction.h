#include <iostream> 
#include <fstream>
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp> 
#include<opencv2/highgui/highgui.hpp> 
#include <opencv2/imgproc/imgproc.hpp>
#include <math.h>
#include <vector>
#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkProperty.h>
#include <vtkTubeFilter.h>
#include <vtkParametricSpline.h>
#include <vtkParametricFunctionSource.h>
#include <vtkSmartPointer.h>
#include <vtkSTLReader.h>
#include <vtkSTLWriter.h>
#include <vtkAppendPolyData.h>
#include <vtkCleanPolyData.h>
#include <vtkBooleanOperationPolyDataFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkTriangleFilter.h>
#include <vtkCylinderSource.h>
#include <vtkTupleInterpolator.h>
#include <vtkDoubleArray.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataMapper.h>

using namespace cv; //ʹ��OpenCV�����ռ�
using namespace std;//ʹ��C++�����ռ�

#define  PI 3.14159265358979323846 //�궨��Բ����PI

typedef struct My_Scalefactor//�����������ӽṹ��
{
	double IamgeAXfactor;//ImageA�϶�ά��X�������D1
	double ImageAYfactor;//ImageA�϶�ά��Y�������D1
	double ImageBXfactor;//ImageB�϶�ά��X�������D2
	double ImageBYfactor;//ImageB�϶�ά��X�������D2
};
typedef struct My2dPoint//�����Լ��Ķ�ά�����ṹ��
{
	double x;
	double y;
};
typedef struct My3dPoint//�����Լ�����ά�����ṹ��
{
	double x;
	double y;
	double z;
};

void InitMat(Mat& m, double* num);  //��һά�����ʼ��OpenCV����Mat����
Mat GenerateRymat(double jiaodu);  //����Ry����OpenCV����Mat����,�����Ƕ�Ϊ�Ƕȣ�ֱ�ӷ���һ��Mat����
Mat GenerateRxmat(double jiaodu);  //����Rx����OpenCV����Mat����,�����Ƕ�Ϊ�Ƕ�,ֱ�ӷ���һ��Mat����
float ExtractMatIJValue(Mat& M, int i, int j);  //��ȡMat�����i�е�j�е���ֵ,MΪMat����ijΪ����ֵ
Mat Generate3DAmat(Mat& M, My_Scalefactor& myfactor);  //�����������ά���A����MΪMat�����,Ϊ�������R����My_ScalefactorΪstruct
Mat Genreate3Damat(Mat& M, My_Scalefactor& myfactor);  //������B�����a����
Mat Generate3DTmat(double T);  //������B�����T����
Mat Generate3Dbmat(Mat& M, My_Scalefactor& myfactor);  //������B�����b����
Mat Generate3DBmat(Mat& a, Mat& b, Mat& t);  //�����������ά���B����
int vtkManual3DReconstruction();  
void vtkManual3DReconstructionRadious();