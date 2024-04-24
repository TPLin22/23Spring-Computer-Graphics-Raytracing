#include "hitsphere.h"
#include "triangle_mesh.h"
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <omp.h>

int main()
{
	std::cout << "������ģ�͵���״��0���֣�1�ǲ��" << std::endl;
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
		// ������ͷ�������߿�����
	float aspect_ratio = 1.5;
	int image_width = 1000;
	int image_height = static_cast<int>(image_width / aspect_ratio);

	// ���ñ�����ɫ
	color background(1.0, 1.0, 1.0);
	//��������
	Scene scene;

	// ����ground���ʺ�λ�ã�ground�Ǹ��뾶1200������,����������
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


	//std::cout << "�������������������ҷ����ڳ����λ���ǹ̶��ģ�������ڴ������޸ġ�\n" << std::endl;
	double r1 = random_double(1.5, 2);
	point3 positionleft = point3(-4, r1, 2);

	//double r2 = 2;
	//point3 positionmiddle = point3(0, r2, -1);

	double r3 = random_double(1.2, 2);
	point3 positionright = point3(4, r3, 1);
	//// ����������ͬ���ʵ�����

	int m1;
	std::cout << "�������һ������Ĳ��ʣ�0�Ǵֲڲ��ʣ�1�ǽ�����2�ǲ���" << std::endl;
	std::cin >> m1;
	if (m1 == 0)
	{
		std::cout << "�ֲڲ���,������" << std::endl;
		// �ֲڱ���,������
		shared_ptr<material> materialrough = make_shared<lambertian>(color::random() * color::random());
		scene.add(make_shared<sphere>(positionleft, r1, materialrough));
	}
	else if (m1 == 1)
	{
		std::cout << "������ȫ����" << std::endl;
		auto materialmetal = make_shared<metal>(color::random(0.3, 1), random_double(0, 0.5));
		scene.add(make_shared<sphere>(positionleft, r1, materialmetal));
	}
	else if (m1 == 2)
	{
		std::cout << "��������" << std::endl;
		shared_ptr<material> materialglass = make_shared<dielectric>(random_double(1, 8.0));
		scene.add(make_shared<sphere>(positionleft, r1, materialglass));
	}
	else
	{
		std::cout << "�ֲڲ���,������" << std::endl;
		shared_ptr<material> materialrough = make_shared<lambertian>(color(0.1, 0.2, 0.9));
		scene.add(make_shared<sphere>(positionleft, r1, materialrough));
	}

	int m2;
	int x2;
	if (m0 == 1) x2 = 0;
	else x2 = -4;
	point3 p2(x2, 0, -1.5);
	std::cout << " ͨ����������ѡ����ʣ�����ֵ������һ������: \n ����Hand/Teapot�Ĳ��ʣ�0:ɢ����ʣ�1:�������ʣ����淴�䣬2:�������ʣ��з�������,Ĭ��ɢ�����" << std::endl;
	std::cin>>m2;
	if(m2==0)
	{
	    std::cout << "�ֲڲ���,������" << std::endl;
	    shared_ptr<material> materialrough = make_shared<lambertian>(color::random() * color::random());
		scene.add(make_shared<triangle_mesh>(filename, materialrough, p2, 0.3));
	}
	else if(m2==1)
	{
	    std::cout << "������ȫ����" << std::endl;
	    auto materialmetal = make_shared<metal>(color::random(0.5, 1), random_double(0, 0.5));
		scene.add(make_shared<triangle_mesh>(filename, materialmetal, p2, 0.3));
	}
	else if(m2==2)
	{
	    std::cout << "��������" << std::endl;
	    shared_ptr<material> materialglass = make_shared<dielectric>(random_double(1, 8.0));
		scene.add(make_shared<triangle_mesh>(filename, materialglass, p2, 0.3));
	}
	else
	{
	    std::cout << "�ֲڲ���,������" << std::endl;
	    shared_ptr<material> materialrough = make_shared<lambertian>(color(0.1, 0.2, 0.9));
		scene.add(make_shared<triangle_mesh>(filename, materialrough, p2, 0.3));
	}

	int m3;
	std::cout << " ͨ����������ѡ����ʣ�����ֵ������һ������: \n ���õڶ�����Ĳ��ʣ�0:ɢ����ʣ�1:�������ʣ����淴�䣬2:�������ʣ��з�������,Ĭ��ɢ�����" << std::endl;
	std::cin >> m3;
	if (m3 == 0)
	{
		std::cout << "�ֲڲ���,������" << std::endl;
		shared_ptr<material> materialrough = make_shared<lambertian>(color::random() * color::random());
		scene.add(make_shared<sphere>(positionright, r3, materialrough));
	}
	else if (m3 == 1)
	{
		std::cout << "������ȫ����" << std::endl;
		auto materialmetal = make_shared<metal>(color::random(0.5, 1), random_double(0, 0.5));
		scene.add(make_shared<sphere>(positionright, r3, materialmetal));
	}
	else if (m3 == 2)
	{
		std::cout << "��������" << std::endl;
		shared_ptr<material> materialglass = make_shared<dielectric>(random_double(1, 8.0));
		scene.add(make_shared<sphere>(positionright, r3, materialglass));
	}
	else
	{
		std::cout << "�ֲڲ���,������" << std::endl;
		shared_ptr<material> materialrough = make_shared<lambertian>(color(0.1, 0.2, 0.9));
		scene.add(make_shared<sphere>(positionright, r3, materialrough));
	}

	int samplesnumber; // ��������
	std::cout << "�������������, Խ�����ʱ��Խ��������ֵ������һ������: ";
	std::cin >> samplesnumber;
	if (samplesnumber == 0)
	{
		samplesnumber = 1;
	}
	int max_depth;// ��������
	std::cout << "��������߷������, Խ�����ʱ��Խ��������ֵ������һ������: ";
	std::cin >> max_depth;
	if (max_depth < 3)
	{
		max_depth = 3;
	}

	// Camera
	point3 lookfrom(0, 3.5, 15);
	point3 lookat(0, 2, 0);
	vec3 vup(0, 1, 0); // view up ������������ϵ����y
	float dist_to_focus = 15.0; //far
	float aperture = 0.25;  // near

	// ��ʼ������Ͳ���
	camera cam(lookfrom, lookat, vup, 42, aspect_ratio, aperture, dist_to_focus);

	std::vector<int>data;
	// �����ͼ��
	cv::Mat result = cv::Mat::zeros(image_height, image_width, CV_8UC3);
	std::cout << "�����У��ȴ�������...." << std::endl;
	double startTime = omp_get_wtime();
	// ʹ��openmp���٣�����12���߳�
	#pragma omp parallel for
	for (int j = image_height - 1; j >= 0; j--)
	{
		for (int i = 0; i < image_width; i++)
		{
			color pixelvale(0, 0, 0);
			for (int s = 0; s < samplesnumber; ++s)
			{
				//��Ӧ�����ӿ�����Ĺ���
				ray r = cam.get_ray(1.0 * i / image_width, 1.0 * j / image_height);
				// ���ػ���
				pixelvale += ray_color(r, scene, background, max_depth);
			}

			result.at<cv::Vec3b>(image_height - j - 1, i)[0] = static_cast<uchar>(255 * clamp(pow(sqrt(1.0 * pixelvale.z() / samplesnumber), 2.0), 0.0, 1.0));
			result.at<cv::Vec3b>(image_height - j - 1, i)[1] = static_cast<uchar>(255 * clamp(pow(sqrt(1.0 * pixelvale.y() / samplesnumber), 2.0), 0.0, 1.0));
			result.at<cv::Vec3b>(image_height - j - 1, i)[2] = static_cast<uchar>(255 * clamp(pow(sqrt(1.0 * pixelvale.x() / samplesnumber), 2.0), 0.0, 1.0));
		}
	}
	double endTime = omp_get_wtime();
	std::cout << "ָ�� 12 ���̣߳�ִ��ʱ��: " << endTime - startTime << " seconds!" << std::endl;
	std::cout << "done !" << std::endl;
	cv::imwrite("resl.png", result);
	std::cout << "����ͼƬ�ɹ�" << std::endl;
	cv::imshow("result", result);
	cv::waitKey(0);

	//int ifloop;
	//std::cout << "�������Ƿ����¼��㣬����0���¼��㣬������������ͼƬ���: ";
	//std::cin >> ifloop;
	//if(ifloop == 0)
	//{
	//    goto LOOP;
	//}
	//else
	//{
	//    cv::imwrite("resl.png", result);
	//    std::cout << "����ͼƬ������" << std::endl;
	//}
}
