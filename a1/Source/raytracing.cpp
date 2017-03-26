// 20.03.17
// code finished, but there are still bugs:
// - some triangle lines still appear
// - try using antialiasing to fix it and extension

#include <iostream>
#include <glm/glm.hpp>
#include <glm/fwd.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"
#include "limits.h"
#include "math.h"

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

vec3 forward(0.0, 0.0, 0.1); //for backward, subtract this vector
vec3 _right(0.1, 0.0, 0.0); //for left, subtract this vector
vec3 up(0.0, -0.1, 0.0); //for down, subtract this vector

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
		cameraPos += forward;
		cout << cameraPos.x << " " << cameraPos.y << " " << cameraPos.z << endl;
	}
	if (keystate[SDLK_DOWN])
	{
		cout << "Backward" << endl;
		//move camera backward
		cameraPos -= forward;
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
		lightPos += forward;
		cout << lightPos.x << " " << lightPos.y << " " << lightPos.z << endl;
	}
	if (keystate[SDLK_s])
	{
		cout << "Light backward" << endl;
		lightPos -= forward;
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

	//cout << cameraPos.x << " " << cameraPos.y << " " << cameraPos.z << " " << endl;

	for( int j=0; j<SCREEN_HEIGHT; ++j )
	{
		for( int i=0; i<SCREEN_WIDTH; ++i )
		{
			vec3 dir((float)i - (float)(SCREEN_WIDTH/2.0),
			       (float)j - (float)(SCREEN_HEIGHT/2.0),
                               (float)focalLength);


			dir = glm::normalize(dir);

			bool check;
			check = ClosestIntersection(cameraPos, dir, triangles, closestInt);

			if (check)
			{
				intersectionIndex = closestInt.triangleIndex;
				//triangle color
				vec3 color(triangles[intersectionIndex].color);
				//PutPixelSDL( screen, i, j, color );

				//light only
				vec3 color2(DirectLight(closestInt));
				//PutPixelSDL( screen, i, j, color2 );

				//color * light
				vec3 color3(color * color2);
				//PutPixelSDL( screen, i, j, color3 );

				//indirect illumination
				vec3 _R = color * (color2 + indirectLight);
				PutPixelSDL(screen, i, j, _R);
			}
			else
			{
				//black
				vec3 color( 0.0, 0.0, 0.0 );
				PutPixelSDL( screen, i, j, color );
			}
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
	//cout << _normal.x << " " << _normal.y << " " << _normal.z << endl;

	//rx = lightPos.x - intersectionPos.x;
	//ry = lightPos.y - intersectionPos.y;
	//rz = lightPos.z - intersectionPos.z;
	//cout << rx << " " << ry << " " << rz << endl;
	//cout << intersectionPos.x << " " << intersectionPos.y << " " << intersectionPos.z << endl;

	//radius = sqrt(pow(rx, 2) + pow(ry, 2) + pow(rz, 2));

	radius = glm::distance(lightPos, intersectionPos);
	area = 4 * PI * radius * radius;
	//cout << radius << " " << area << endl;

	//_n = glm::normalize(_n);
	_r = lightPos - surfacePoint;
	_r = glm::normalize(_r);
	//cout << _r.x << " " << _r.y << " " << _r.z << endl;

	float pMax = glm::max((float)glm::dot(_r, _n), (float)0.0);
	//cout << pMax << endl;

	//d = lightColor/area;

	if (pMax > 0.0)
	{
		d = lightColor/area;
	}
	else
	{
		d = vec3(0.0, 0.0, 0.0);
	}

	//cout << d.x << " " << d.y << " " << d.z << endl;

	/* check for another surface */
	/*
	bool check = ClosestIntersection(intersectionPos, lightPos, triangles, otherSurface);
	if (check)
	{
		if (glm::distance(otherSurface.position, intersectionPos) < radius)
		{
			d = vec3(0.0, 0.0, 0.0);
		}
	}
	*/

	return d;
}

float findMax(float a, float b)
{
	if (a > b)
		return a;
	else
		return b;
}
