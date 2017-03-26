//assignment 0.2

#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"

using namespace std;
using glm::vec3;
using glm::vec2;
using glm::mat3;

/* ----------------------------------------------------------------------------*/
/* GLOBAL VARIABLES                                                            */

const int SCREEN_WIDTH = 500;
const int SCREEN_HEIGHT = 500;
SDL_Surface* screen;
int t;
vector<vec3> stars(1000);
vector<vec2> stars2D(1000);

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

void Update();
void Draw();

int main( int argc, char* argv[] )
{
	screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT );
	t = SDL_GetTicks();	// Set start value for timer.
	
	for (size_t i=0; i<stars.size(); ++i)
	{
		float rx = float(rand())/float(RAND_MAX);
		float ry = float(rand())/float(RAND_MAX);
		float rz = float(rand())/float(RAND_MAX);

		//min + r * (max - min)
		float xr = -1 + rx * 2;
		float yr = -1 + ry * 2;
		float zr = rz;
		
		stars[i].x = xr;
		stars[i].y = yr;
		stars[i].z = zr;
	}

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
	
	float vel = 0.001;

	for (size_t s=0; s<stars.size(); ++s)
	{
		stars[s].z = stars[s].z - (vel * dt);
		
		if (stars[s].z <= 0)
			stars[s].z += 1;
		if (stars[s].z > 1)
			stars[s].z -= 1;
	}
}

void Draw()
{
	SDL_FillRect(screen, 0, 0);

	if( SDL_MUSTLOCK(screen) )
		SDL_LockSurface(screen);

	float foc = SCREEN_HEIGHT/2;
	vec3 color(1,1,1);
	
	for (size_t s=0; s<stars.size(); ++s)
	{
		stars2D[s].x = (foc*(stars[s].x/stars[s].z)) + SCREEN_WIDTH/2;
		stars2D[s].y = (foc*(stars[s].y/stars[s].z)) + SCREEN_HEIGHT/2;

		vec3 color2 = 0.2f * vec3(1,1,1)/(stars[s].z * stars[s].z);
		PutPixelSDL( screen, stars2D[s].x, stars2D[s].y, color2 );
	}

	if( SDL_MUSTLOCK(screen) )
		SDL_UnlockSurface(screen);

	SDL_UpdateRect( screen, 0, 0, 0, 0 );
}
