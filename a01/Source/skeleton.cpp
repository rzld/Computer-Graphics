#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"
#include <vector>

using namespace std;
using glm::vec3;
using glm::mat3;

/* ----------------------------------------------------------------------------*/
/* GLOBAL VARIABLES                                                            */

const int SCREEN_WIDTH = 500;
const int SCREEN_HEIGHT = 500;
SDL_Surface* screen;
int t;

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

void Update();
void Draw();
void Interpolate(float a, float b, vector<float>& result);
void Interpolate2(vec3 a, vec3 b, vector<vec3>& result);

int main( int argc, char* argv[] )
{
	screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT );
	t = SDL_GetTicks();	// Set start value for timer.
	
	//vector<float> result(5);
	//Interpolate(5, 14, result);
	//for (int i=0; i<int(result.size()); ++i)
	//	cout << result[i] << " ";

//	vector<vec3> result(4);
//	vec3 a(1, 4, 9.2);
//	vec3 b(4, 1, 9.8);
//	Interpolate2(a, b, result);
//	for (int i=0; i<int(result.size()); ++i)
//	{
//		cout << "( " 
//		<< result[i].x << ", " 
//		<< result[i].y << ", " 
//		<< result[i].z << " ) ";
//	}

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
	//cout << "Render time: " << dt << " ms." << endl;
}

void Draw()
{
	if( SDL_MUSTLOCK(screen) )
		SDL_LockSurface(screen);
	
	vec3 topLeft(1,0,1); //red
	vec3 topRight(0,0,1); //blue
	vec3 bottomRight(0,1,0); //green
	vec3 bottomLeft(1,1,0); //yellow

	vector<vec3> leftSide(SCREEN_HEIGHT);
	vector<vec3> rightSide(SCREEN_HEIGHT);
	Interpolate2(topLeft, bottomLeft, leftSide);
	Interpolate2(topRight, bottomRight, rightSide);
	
	vector<vec3> interpolatedRow(SCREEN_WIDTH);

	for( int b=0; b<SCREEN_HEIGHT; ++b )
	{
		Interpolate2(leftSide[b], rightSide[b], interpolatedRow);
		for( int a=0; a<SCREEN_WIDTH; ++a )
		{
			vec3 color( 0.0, 1.0, 0.5 );
			vec3 color2( interpolatedRow[a].x, interpolatedRow[a].y, interpolatedRow[a].z );
			PutPixelSDL( screen, a, b, color2 );
		}
	}

	if( SDL_MUSTLOCK(screen) )
		SDL_UnlockSurface(screen);

	SDL_UpdateRect( screen, 0, 0, 0, 0 );
}

void Interpolate(float a, float b, vector<float>& result)
{
	int vecSize = result.size();	
	float d = b - a;
	float e = d / (vecSize-1);
	result[0] = a;
	for (int i = 1; i < vecSize; i++)
	{
		float newValue = result[i-1] + e;
		result[i] = newValue;
	}
}

void Interpolate2(vec3 a, vec3 b, vector<vec3>& result)
{
	int vecSize = result.size();
	vector<float> xResult(vecSize), yResult(vecSize), zResult(vecSize);
	
	Interpolate(a.x, b.x, xResult);
	Interpolate(a.y, b.y, yResult);
	Interpolate(a.z, b.z, zResult);

	result[0].x = a.x;
	result[0].y = a.y;
	result[0].z = a.z;

	for (int i = 1; i < vecSize; i++)
	{
		result[i].x = xResult[i];
		result[i].y = yResult[i];
		result[i].z = zResult[i];
	}
}

