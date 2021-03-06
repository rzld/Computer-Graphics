//this file is for raytracer with extension
// extension: anti-aliasing (smooth edges)

#include <iostream>
#include <glm/glm.hpp>
#include <glm/fwd.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"
#include "limits.h"
#include "math.h"
//#include <random>
#include <fstream>

#define PI 3.14159

using namespace std;
using glm::vec3;
using glm::mat3;

/* ----------------------------------------------------------------------------*/
/* GLOBAL VARIABLES                                                            */

struct Intersection
{
	vec3 position;
	float distance;
	int triangleIndex;
};

const int SCREEN_WIDTH = 500;
const int SCREEN_HEIGHT = 500;
SDL_Surface* screen;
int t;

vector<Triangle> triangles;
Intersection closestInt;
vec3 s;

float focalLength = 250;
vec3 cameraPos(0.0, 0.0, -1.9);

float m = std::numeric_limits<float>::max();
int intersectionIndex;

mat3 R; //camera rotation matrix
float yaw; //y-axis angle for rotation

vec3 lightPos(0.0, -0.5, -0.5);
vec3 lightColor = 14.f * vec3(1.0, 1.0, 1.0);
vec3 indirectLight = 0.5f * vec3(1, 1, 1);

vec3 forwards(0.0, 0.0, 0.1); //for backward, subtract this vector
vec3 _right(0.1, 0.0, 0.0); //for left, subtract this vector
vec3 up(0.0, -0.1, 0.0); //for down, subtract this vector

//DOF
int apertureSize = 4;
const int rays = 8;
float focalDistance = 300;

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

void Update();
void Draw();
bool ClosestIntersection(vec3 start, vec3 dir, const vector<Triangle>& triangles, Intersection& closestIntersection);
vec3 DirectLight(const Intersection& i);
float findMax(float a, float b);

int main( int argc, char* argv[] )
{
	screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT );
	t = SDL_GetTicks();	// Set start value for timer.

	LoadTestModel(triangles);
	cout << "Load test model" << endl;
	//cout << lightColor.x << " " << lightColor.y << " " << lightColor.z << endl;

	while( NoQuitMessageSDL() )
	{
		Update();
		Draw();
	}

	SDL_SaveBMP( screen, "screenshot.bmp" );
	return 0;
}

void Update()
{
	// Compute frame time:
	int t2 = SDL_GetTicks();
	float dt = float(t2-t);
	t = t2;
	cout << "Render time: " << dt << " ms." << endl;

	//CAMERA CONTROL

	//vec3 right(	R[0][0], R[0][1], R[0][2]	);
	//vec3 down(	R[1][0], R[1][1], R[1][2]	);
	//vec3 forward(	R[2][0], R[2][1], R[2][2]	);

	Uint8* keystate = SDL_GetKeyState(0);
	if (keystate[SDLK_UP])
	{
		cout << "Forward" << endl;
		//move camera forward
		cameraPos += forwards;
		cout << cameraPos.x << " " << cameraPos.y << " " << cameraPos.z << endl;
	}
	if (keystate[SDLK_DOWN])
	{
		cout << "Backward" << endl;
		//move camera backward
		cameraPos -= forwards;
		cout << cameraPos.x << " " << cameraPos.y << " " << cameraPos.z << endl;
	}
	if (keystate[SDLK_LEFT])
	{
		cout << "Left" << endl;
		//move camera to the left
		cameraPos -= _right;
		cout << cameraPos.x << " " << cameraPos.y << " " << cameraPos.z << endl;
		//rotate camera to the left
		//yaw = ;

	}
	if (keystate[SDLK_RIGHT])
	{
		cout << "Right" << endl;
		//move camera to the right
		cameraPos += _right;
		cout << cameraPos.x << " " << cameraPos.y << " " << cameraPos.z << endl;
		//rotate camera to the right
		//yaw = ;
	}
	if (keystate[SDLK_o])
	{
		cout << "Up" << endl;
		//move camera to the right
		cameraPos += up;
		cout << cameraPos.x << " " << cameraPos.y << " " << cameraPos.z << endl;
		//rotate camera to the right
		//yaw = ;
	}
	if (keystate[SDLK_l])
	{
		cout << "Down" << endl;
		//move camera to the right
		cameraPos -= up;
		cout << cameraPos.x << " " << cameraPos.y << " " << cameraPos.z << endl;
		//rotate camera to the right
		//yaw = ;
	}


	//LIGHT CONTROL

	if (keystate[SDLK_w])
	{
		cout << "Light forward" << endl;
		lightPos += forwards;
		cout << lightPos.x << " " << lightPos.y << " " << lightPos.z << endl;
	}
	if (keystate[SDLK_s])
	{
		cout << "Light backward" << endl;
		lightPos -= forwards;
		cout << lightPos.x << " " << lightPos.y << " " << lightPos.z << endl;
	}
	if (keystate[SDLK_a])
	{
		cout << "Light left" << endl;
		lightPos -= _right;
		cout << lightPos.x << " " << lightPos.y << " " << lightPos.z << endl;
	}
	if (keystate[SDLK_d])
	{
		cout << "Light right" << endl;
		lightPos += _right;
		cout << lightPos.x << " " << lightPos.y << " " << lightPos.z << endl;
	}
	if (keystate[SDLK_q])
	{
		cout << "Light up" << endl;
		lightPos += up;
		cout << lightPos.x << " " << lightPos.y << " " << lightPos.z << endl;
	}
	if (keystate[SDLK_e])
	{
		cout << "Light down" << endl;
		lightPos -= up;
		cout << lightPos.x << " " << lightPos.y << " " << lightPos.z << endl;
	}

}

void Draw()
{
	if( SDL_MUSTLOCK(screen) )
		SDL_LockSurface(screen);

	float jitterMat[8] = { -1.0 / 4.0,  3.0 / 4.0,
						   3.0 / 4.0,  1.0 / 4.0,
						  -3.0 / 4.0, -1.0 / 4.0,
						   1.0 / 4.0, -3.0 / 4.0 };

	//cout << cameraPos.x << " " << cameraPos.y << " " << cameraPos.z << " " << endl;
	float eyeDist = glm::distance(cameraPos, vec3(0.0, 0.0, float(focalLength)));
	for( int j=0; j<SCREEN_HEIGHT; ++j )
	{
		for( int i=0; i<SCREEN_WIDTH; ++i )
		{
			vector<vec3> newCameraPos(rays), newCameraPosN(rays);

			vec3 dir((float)i - (float)(SCREEN_WIDTH/2.0),
			         (float)j - (float)(SCREEN_HEIGHT/2.0),
                     (float)focalLength);

			vec3 dir2 = glm::normalize(dir);

			float eyeToPixel = glm::distance(cameraPos, dir);

			vec3 focalPoint = cameraPos + (eyeToPixel/(eyeDist/(eyeDist+focalDistance))) * (dir-cameraPos);

			float randomX[rays], randomY[rays], _x, _y;
			//std::default_random_engine generator;
  			//std::uniform_int_distribution<int> distribution(0,apertureSize);

			//for (int n=0; n<rays; n++)
			//{
			//	//_x = distribution(generator);
			//	//_y = distribution(generator);

			//	_x = 0 + (rand()/(RAND_MAX/(apertureSize-0)));
			//	_y = 0 + (rand()/(RAND_MAX/(apertureSize-0)));
			//	randomX[n] = _x;
			//	randomY[n] = _y;

			//	//new camera position: in random, with aperture size
			//	newCameraPos[n] = vec3(dir.x + (float)randomX[n], 
			//						   dir.y + (float)randomY[n], 
			//							(float)focalLength);
			//	newCameraPosN[n] = glm::normalize(newCameraPos[n]);

			//	/*newCameraPos[n] = vec3((dir.x - apertureSize / 2) + n,
			//						   (dir.y - apertureSize / 2) + n,
			//						   (float)focalLength);

			//	newCameraPosN[n] = glm::normalize(newCameraPos[n]);*/
			//	//myfile << newCameraPos[n].x << " " << newCameraPos[n].y << " " << newCameraPos[n].z << endl;
			//}

			for (int n = 0; n < 4; ++n)
			{
				newCameraPos[n] = vec3(dir.x + jitterMat[2 * n],
					dir.y + jitterMat[2 * n + 1],
					dir.z);
				
				newCameraPosN[n] = glm::normalize(newCameraPos[n]);
			}

			vec3 newDir = glm::normalize(focalPoint);
			
			vec3 totalColor;
			int nx = 0;

			for (int n=0; n<4; n++)
			{
				bool check;
				check = ClosestIntersection(cameraPos, newCameraPosN[n], triangles, closestInt);
				vec3 color;

				if (check)
				{
					//cout << "Check" << endl;
					intersectionIndex = closestInt.triangleIndex;
					//triangle color
					vec3 colorOri = vec3(triangles[intersectionIndex].color);
					//PutPixelSDL( screen, i, j, color );

					//light only
					vec3 lightVec(DirectLight(closestInt));
					//PutPixelSDL( screen, i, j, color2 );

					//color * light
					//vec3 colorLight(colorOri * lightVec);
					//PutPixelSDL( screen, i, j, color3 );

					//indirect illumination
					color = colorOri * (lightVec + indirectLight);
					//PutPixelSDL(screen, i, j, _R);
				}
				else
				{
					//black
					color = vec3( 0.0, 0.0, 0.0 );
					//PutPixelSDL( screen, i, j, color );
				}

				if (color.x > 0.0 && color.y > 0.0 && color.z > 0.0)
				{
					totalColor += color;
					nx++;
				}
			}
			//cout << totalColor.x << " " << totalColor.y << " " << totalColor.z << endl;
			totalColor /= (float)nx;

			PutPixelSDL(screen, i, j, totalColor);
		}
	}

	if( SDL_MUSTLOCK(screen) )
		SDL_UnlockSurface(screen);

	SDL_UpdateRect( screen, 0, 0, 0, 0 );
}

bool ClosestIntersection(vec3 start, vec3 dir, const vector<Triangle>& triangles, Intersection& closestIntersection)
{
	float _t, _u, _v;

	closestIntersection.distance = 99.0;
	int tSize = triangles.size();
	//float dist[tSize];
	//td::vector<vec3> pos;
	//cout << tSize << endl;

	for (int i = 0; i < tSize; ++i)
	{
		//cout << "Triangle #" << i << endl;

		vec3 v0 = triangles[i].v0;
		vec3 v1 = triangles[i].v1;
		vec3 v2 = triangles[i].v2;

		//cout << "checkpoint 1 #" << i << endl;
		vec3 e1 = v1 - v0;
		vec3 e2 = v2 - v0;
		vec3 b = start - v0;

		//cout << "checkpoint 3 #" << i << endl;
		mat3 A(-dir, e1, e2);
		vec3 xx = glm::inverse(A)*b;

		_t = xx.x;
		_u = xx.y;
		_v = xx.z;
		//dist[i] = _t;

		if (_u > 0.0 && _v > 0.0 && _u + _v < 1.0 && _t < closestIntersection.distance && _t >= 0.0)
		{
				//cout << "check " << dist[i] << endl;
				closestIntersection.distance = _t;
				closestIntersection.position = start + _t * dir;
				closestIntersection.triangleIndex = i;
		}
	}

	//intersectionIndex = closestIntersection.triangleIndex;
	//cout << intersectionIndex << endl;
	if (closestIntersection.distance == 99.0)
		return false;
	else
		return true;

}

vec3 DirectLight(const Intersection& i)
{
	float radius;
	float area;
	vec3 _n, _r, d;
	vec3 intersectionPos(i.position);
	vec3 surfacePoint(triangles[i.triangleIndex].v0);
	Intersection otherSurface;

	_n = triangles[i.triangleIndex].normal;

	radius = glm::distance(lightPos, intersectionPos);
	area = (float)(4.0 * PI * radius * radius);

	_r = lightPos - intersectionPos;
	_r = glm::normalize(_r);

	float pMax = glm::max((float)glm::dot(_r, _n), (float)0.0);

	if (pMax > 0.0)
	{
		d = lightColor/area;
	}
	else
	{
		d = vec3(0.0, 0.0, 0.0);
	}

	/* check for another surface */

	bool check = ClosestIntersection(intersectionPos*0.99f, _r, triangles, otherSurface);
	if (check && glm::distance(intersectionPos*0.99f, otherSurface.position) < radius)
	{
		//cout << "Shadow!" << endl;
		d = vec3(0.0, 0.0, 0.0);
	}


	return d;
}

float findMax(float a, float b)
{
	if (a > b)
		return a;
	else
		return b;
}
