#include "hitsphere.h"
#include "triangle_mesh.h"
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <omp.h>

int main()
{
	std::cout << "请输入模型的形状，0是手，1是茶壶" << std::endl;
	int m0;
	std::cin >> m0;
	//char filename[] = "teapot.obj";
	//char filename[] = "Hand.fbx";
	char filename[100];

	if (m0 == 1)
	{
		strcpy_s(filename, "teapot.obj");
	}
	else
	{
		strcpy_s(filename, "Hand.fbx");
	}
	
	//LOOP:
		// 设摄像头画幅，高宽，比例
	float aspect_ratio = 1.5;
	int image_width = 1000;
	int image_height = static_cast<int>(image_width / aspect_ratio);

	// 设置背景颜色
	color background(1.0, 1.0, 1.0);
	//创建场景
	Scene scene;

	// 设置ground材质和位置，ground是个半径1200的球体,表面漫反射
	shared_ptr<material> groundmaterial = make_shared<lambertian>(color(0.6, 0.5, 0.52));
	shared_ptr<material> wallmaterial = make_shared<lambertian>(color(1.02 / 2.55, 2.04 / 2.55, 1));
	shared_ptr<material> wallmaterial2 = make_shared<lambertian>(color(1.48 / 2.55, 0.07 / 2.55, 0.1 / 2.55));
	shared_ptr<material> wallmaterial3 = make_shared<lambertian>(color(0.67 / 2.55, 1.28 / 2.55, 1.09 / 2.55));
	shared_ptr<material> topmaterial = make_shared<glittering>(color(1, 1, 1));
	shared_ptr<material> metalmaterial = make_shared<metal>(color(0.8, 0.8, 0.8), 0.3);

	scene.add(make_shared<plane>(point3(2, 0, 2), vec3(0, 1, 0), groundmaterial));
	scene.add(make_shared<plane>(point3(8, 0, 0), vec3(-1, 0, 0), wallmaterial));
	scene.add(make_shared<plane>(point3(-8, 0, 0), vec3(1, 0, 0), wallmaterial3));
	scene.add(make_shared<plane>(point3(-10, 0, -10), vec3(0, 0, 1), wallmaterial2));
	scene.add(make_shared<plane>(point3(2, 8, 2), vec3(0, -1, 0), topmaterial));


	//std::cout << "场景中有三个球，左中右放置在场景里，位置是固定的，你可以在代码里修改。\n" << std::endl;
	double r1 = random_double(1.5, 2);
	point3 positionleft = point3(-4, r1, 2);

	//double r2 = 2;
	//point3 positionmiddle = point3(0, r2, -1);

	double r3 = random_double(1.2, 2);
	point3 positionright = point3(4, r3, 1);
	//// 设置三个不同材质的球体

	int m1;
	std::cout << "请输入第一个球体的材质，0是粗糙材质，1是金属，2是玻璃" << std::endl;
	std::cin >> m1;
	if (m1 == 0)
	{
		std::cout << "粗糙材质,漫反射" << std::endl;
		// 粗糙表面,漫反射
		shared_ptr<material> materialrough = make_shared<lambertian>(color::random() * color::random());
		scene.add(make_shared<sphere>(positionleft, r1, materialrough));
	}
	else if (m1 == 1)
	{
		std::cout << "金属，全反射" << std::endl;
		auto materialmetal = make_shared<metal>(color::random(0.3, 1), random_double(0, 0.5));
		scene.add(make_shared<sphere>(positionleft, r1, materialmetal));
	}
	else if (m1 == 2)
	{
		std::cout << "玻璃材质" << std::endl;
		shared_ptr<material> materialglass = make_shared<dielectric>(random_double(1, 8.0));
		scene.add(make_shared<sphere>(positionleft, r1, materialglass));
	}
	else
	{
		std::cout << "粗糙材质,漫反射" << std::endl;
		shared_ptr<material> materialrough = make_shared<lambertian>(color(0.1, 0.2, 0.9));
		scene.add(make_shared<sphere>(positionleft, r1, materialrough));
	}

	int m2;
	int x2;
	if (m0 == 1) x2 = 0;
	else x2 = -4;
	point3 p2(x2, 0, -1.5);
	std::cout << " 通过输入数字选择材质，输入值必须是一个数字: \n 设置Hand/Teapot的材质，0:散射材质，1:金属材质，表面反射，2:玻璃材质，有反射折射,默认散射材质" << std::endl;
	std::cin>>m2;
	if(m2==0)
	{
	    std::cout << "粗糙材质,漫反射" << std::endl;
	    shared_ptr<material> materialrough = make_shared<lambertian>(color::random() * color::random());
		scene.add(make_shared<triangle_mesh>(filename, materialrough, p2, 0.3));
	}
	else if(m2==1)
	{
	    std::cout << "金属，全反射" << std::endl;
	    auto materialmetal = make_shared<metal>(color::random(0.5, 1), random_double(0, 0.5));
		scene.add(make_shared<triangle_mesh>(filename, materialmetal, p2, 0.3));
	}
	else if(m2==2)
	{
	    std::cout << "玻璃材质" << std::endl;
	    shared_ptr<material> materialglass = make_shared<dielectric>(random_double(1, 8.0));
		scene.add(make_shared<triangle_mesh>(filename, materialglass, p2, 0.3));
	}
	else
	{
	    std::cout << "粗糙材质,漫反射" << std::endl;
	    shared_ptr<material> materialrough = make_shared<lambertian>(color(0.1, 0.2, 0.9));
		scene.add(make_shared<triangle_mesh>(filename, materialrough, p2, 0.3));
	}

	int m3;
	std::cout << " 通过输入数字选择材质，输入值必须是一个数字: \n 设置第二个球的材质，0:散射材质，1:金属材质，表面反射，2:玻璃材质，有反射折射,默认散射材质" << std::endl;
	std::cin >> m3;
	if (m3 == 0)
	{
		std::cout << "粗糙材质,漫反射" << std::endl;
		shared_ptr<material> materialrough = make_shared<lambertian>(color::random() * color::random());
		scene.add(make_shared<sphere>(positionright, r3, materialrough));
	}
	else if (m3 == 1)
	{
		std::cout << "金属，全反射" << std::endl;
		auto materialmetal = make_shared<metal>(color::random(0.5, 1), random_double(0, 0.5));
		scene.add(make_shared<sphere>(positionright, r3, materialmetal));
	}
	else if (m3 == 2)
	{
		std::cout << "玻璃材质" << std::endl;
		shared_ptr<material> materialglass = make_shared<dielectric>(random_double(1, 8.0));
		scene.add(make_shared<sphere>(positionright, r3, materialglass));
	}
	else
	{
		std::cout << "粗糙材质,漫反射" << std::endl;
		shared_ptr<material> materialrough = make_shared<lambertian>(color(0.1, 0.2, 0.9));
		scene.add(make_shared<sphere>(positionright, r3, materialrough));
	}

	int samplesnumber; // 采样次数
	std::cout << "请输入采样次数, 越大计算时间越长，输入值必须是一个数字: ";
	std::cin >> samplesnumber;
	if (samplesnumber == 0)
	{
		samplesnumber = 1;
	}
	int max_depth;// 最大反射次数
	std::cout << "请输入光线反射次数, 越大计算时间越长，输入值必须是一个数字: ";
	std::cin >> max_depth;
	if (max_depth < 3)
	{
		max_depth = 3;
	}

	// Camera
	point3 lookfrom(0, 3.5, 15);
	point3 lookat(0, 2, 0);
	vec3 vup(0, 1, 0); // view up 方向，右手坐标系，正y
	float dist_to_focus = 15.0; //far
	float aperture = 0.25;  // near

	// 初始化相机和参数
	camera cam(lookfrom, lookat, vup, 42, aspect_ratio, aperture, dist_to_focus);

	std::vector<int>data;
	// 输出的图像
	cv::Mat result = cv::Mat::zeros(image_height, image_width, CV_8UC3);
	std::cout << "计算中，等待结果输出...." << std::endl;
	double startTime = omp_get_wtime();
	// 使用openmp加速，设置12个线程
	#pragma omp parallel for
	for (int j = image_height - 1; j >= 0; j--)
	{
		for (int i = 0; i < image_width; i++)
		{
			color pixelvale(0, 0, 0);
			for (int s = 0; s < samplesnumber; ++s)
			{
				//对应像素视口射出的光线
				ray r = cam.get_ray(1.0 * i / image_width, 1.0 * j / image_height);
				// 像素积分
				pixelvale += ray_color(r, scene, background, max_depth);
			}

			result.at<cv::Vec3b>(image_height - j - 1, i)[0] = static_cast<uchar>(255 * clamp(pow(sqrt(1.0 * pixelvale.z() / samplesnumber), 2.0), 0.0, 1.0));
			result.at<cv::Vec3b>(image_height - j - 1, i)[1] = static_cast<uchar>(255 * clamp(pow(sqrt(1.0 * pixelvale.y() / samplesnumber), 2.0), 0.0, 1.0));
			result.at<cv::Vec3b>(image_height - j - 1, i)[2] = static_cast<uchar>(255 * clamp(pow(sqrt(1.0 * pixelvale.x() / samplesnumber), 2.0), 0.0, 1.0));
		}
	}
	double endTime = omp_get_wtime();
	std::cout << "指定 12 个线程，执行时间: " << endTime - startTime << " seconds!" << std::endl;
	std::cout << "done !" << std::endl;
	cv::imwrite("resl.png", result);
	std::cout << "保存图片成功" << std::endl;
	cv::imshow("result", result);
	cv::waitKey(0);

	//int ifloop;
	//std::cout << "请输入是否重新计算，输入0重新计算，输入其他保存图片结果: ";
	//std::cin >> ifloop;
	//if(ifloop == 0)
	//{
	//    goto LOOP;
	//}
	//else
	//{
	//    cv::imwrite("resl.png", result);
	//    std::cout << "保存图片并结束" << std::endl;
	//}
}
