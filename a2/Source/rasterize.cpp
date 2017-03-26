//RASTERISATION

//Progress 11/03/2017:
// Still confused about zinv

#include <iostream>
#include <glm/glm.hpp>
#include <glm/fwd.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"
#include <limits>

using namespace std;
using glm::vec3;
using glm::mat3;
using glm::vec2;

/* ----------------------------------------------------------------------------*/
/* GLOBAL VARIABLES                                                            */

const int SCREEN_WIDTH = 500;
const int SCREEN_HEIGHT = 500;
SDL_Surface* screen;
int t;

vector<Triangle> triangles;
vec3 cameraPos(0, 0, -3.001);
float focalLength = 250.0;
mat3 R;
float yaw = 0;
float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];
vec3 currentColor;
vec3 lightPos(0, -0.5, -0.7);
vec3 lightPower = 1.1f*vec3(1,1,1);
vec3 indirectLight = 0.5f*vec3(1,1,1);

struct Pixel
{
	int x;
	int y;
	float zinv;
};
struct Vertex
{
	vec3 position;
	vec3 normal;
	vec3 reflectance;
};

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

void Update();
void Draw();
void VertexShader(const vec3& v, glm::ivec2& p);
void Interpolate(glm::ivec2 a, glm::ivec2 b, vector<glm::ivec2>& result);
void DrawLineSDL(SDL_Surface* surface, glm::ivec2 a, glm::ivec2 b, vec3 color);
void DrawPolygonEdges(const vector<vec3>& vertices);
void ComputePolygonRows(const vector<glm::ivec2>& vertexPixels,
			      vector<glm::ivec2>& leftPixels,
			      vector<glm::ivec2>& rightPixels);
void DrawPolygonRows(const vector<glm::ivec2>& leftPixels, const vector<glm::ivec2>& rightPixels);
void DrawPolygon(const vector<vec3>& vertices);

void InterpolateP(Pixel a, Pixel b, vector<Pixel>& result);
void ComputePolygonRowsP(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels);
void DrawPolygonRowsP(const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels);
void VertexShaderP(const vec3& v, Pixel& p);
void PixelShader(const Pixel& p);

int main( int argc, char* argv[] )
{
	screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT );
	t = SDL_GetTicks();	// Set start value for timer.

	LoadTestModel(triangles);
	cout << "Load test model" << endl;

	// ComputePolygonRows test
	vector<glm::ivec2> vertexPixels(3);
	vertexPixels[0] = glm::ivec2(10, 5);
	vertexPixels[1] = glm::ivec2(5, 10);
	vertexPixels[2] = glm::ivec2(15, 15);
	vector<glm::ivec2> leftPixels;
	vector<glm::ivec2> rightPixels;
	ComputePolygonRows(vertexPixels, leftPixels, rightPixels);
	for (int row=0; row<leftPixels.size(); ++row)
	{
		cout << "Start: ("
		     << leftPixels[row].x << ","
		     << leftPixels[row].y << "). "
		     << "End: ("
		     << rightPixels[row].x << ","
		     << rightPixels[row].y << "). " << endl;
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
	//cout << "Render time: " << dt << " ms." << endl;
}

void Draw()
{
	SDL_FillRect(screen,0,0);

	if( SDL_MUSTLOCK(screen) )
		SDL_LockSurface(screen);

	for (int y=0; y<SCREEN_HEIGHT; ++y)
		for (int x=0; x<SCREEN_WIDTH; ++x)
			depthBuffer[y][x] = 0;

	for (int i=0; i<triangles.size(); ++i)
	{
		vector<vec3> vertices(3);

		vertices[0] = triangles[i].v0;
		vertices[1] = triangles[i].v1;
		vertices[2] = triangles[i].v2;

		// dots only
		for (int v=0; v<3; ++v)
		{
			glm::ivec2 projPos;
			VertexShader(vertices[v], projPos);
			//cout << projPos.x << " " << projPos.y << endl;
			//vec3 color(1.0, 1.0, 1.0);
			//PutPixelSDL(screen, projPos.x, projPos.y, color);
		}

		DrawPolygonEdges(vertices);
		//DrawPolygon(vertices);
	}

	//original
	/*for( int y=0; y<SCREEN_HEIGHT; ++y )
	{
		for( int x=0; x<SCREEN_WIDTH; ++x )
		{
			vec3 color( 1.0, 0.0, 0.0 );
			PutPixelSDL( screen, x, y, color );
		}
	}*/

	if( SDL_MUSTLOCK(screen) )
		SDL_UnlockSurface(screen);

	SDL_UpdateRect( screen, 0, 0, 0, 0 );
}
void VertexShader(const vec3& v, glm::ivec2& p)
{
	// take 3D position of a vertex v and compute its 2D img position and store it in p
	vec3 pos = (v - cameraPos) * R;

	p.x = (focalLength * ( pos.x /  pos.z )) + (SCREEN_WIDTH / 2);
	p.y = (focalLength * ( pos.y /  pos.z )) + (SCREEN_HEIGHT / 2);
}

void Interpolate(glm::ivec2 a, glm::ivec2 b, vector<glm::ivec2>& result)
{
	int N = result.size();
	vec2 step = vec2(b-a) / float(max(N-1,1));
	vec2 current(a);

	for (int i=0; i<N; ++i)
	{
		result[i] = current;
		current += step;
	}
}

//draw line
void DrawLineSDL(SDL_Surface* surface, glm::ivec2 a, glm::ivec2 b, vec3 color)
{
	glm::ivec2 delta = glm::abs(a-b);
	int pixels = glm::max(delta.x, delta.y) + 1;

	vector<glm::ivec2> line(pixels);
	Interpolate(a, b, line);

	for (int i=0; i<line.size(); i++)
	{
		PutPixelSDL(surface, line[i].x, line[i].y, color);
	}
}

//draw edges of triangle
void DrawPolygonEdges(const vector<vec3>& vertices)
{
	int V = vertices.size();

	//Transform each vertex from 3D world position to 2D image position
	vector<glm::ivec2> projectedVertices(V);
	for (int i=0; i<V; i++)
	{
		VertexShader(vertices[i], projectedVertices[i]);
	}

	//Loop over all vertices and draw the edge from it to the next vertex
	for (int i=0; i<V; i++)
	{
		int j = (i+1)%V; //next vertex
		vec3 color(1, 1, 1);
		DrawLineSDL(screen, projectedVertices[i], projectedVertices[j], color);
	}
}

//input: projPos
void ComputePolygonRows(const vector<glm::ivec2>& vertexPixels,
			      vector<glm::ivec2>& leftPixels,
			      vector<glm::ivec2>& rightPixels)
{
	int minY = 999;
	int maxY = -999;
	int rows;

	// 1. Find max & min y-value of the polygon and compute the number of rows it occupies
	for (int i=0; i<vertexPixels.size(); i++)
	{
		minY = glm::min(minY, vertexPixels[i].y);
		maxY = glm::max(maxY, vertexPixels[i].x);
	}
	rows = maxY - minY + 1;

	// 2. Resize leftPixels and rightPixels so that they have an element for each row
	leftPixels.resize(rows);
	rightPixels.resize(rows);

	// 3. Initialize the x-coordinates in leftPixels to some really large value
	//    and the x-coordinates in the rightPixels to some really small values
	for (int i=0; i<rows; i++)
	{
		leftPixels[i].x = +numeric_limits<int>::max();
		rightPixels[i].x = -numeric_limits<int>::max();

		leftPixels[i].y = minY + i;
		rightPixels[i].y = minY + i;
	}

	// 4. Loop through all edges of the polygon and use linear interpolation to find
	//    the x-coordinate for each row it occupies. Update the corresponding values in rightPixels
	//    and leftPixels
	for (int i=0; i<vertexPixels.size(); i++)
	{
		int j = (i+1)%vertexPixels.size(); //next vertex
		glm::ivec2 delta = glm::abs(vertexPixels[i]-vertexPixels[j]);
		int pixels = glm::max(delta.x, delta.y) + 1;

		vector<glm::ivec2> line(pixels);
		Interpolate(vertexPixels[i], vertexPixels[j], line);

		// check for every point in line
		for (int a=0; a<line.size(); a++)
		{
			// check for every point in left & right pixels
			for (int b=0; b<rows; b++)
			{
				// if they are in the same y line
				if (leftPixels[b].y == line[a].y)
				{
					int leftX = glm::min(leftPixels[b].x, line[a].x);
					int rightX = glm::max(rightPixels[b].x, line[a].x);

					leftPixels[b].x = leftX;
					rightPixels[b].x = rightX;
				}
			}
		}
	}
}

void DrawPolygonRows(const vector<glm::ivec2>& leftPixels, const vector<glm::ivec2>& rightPixels)
{
	vector<glm::ivec2> rowPixels((rightPixels.x - leftPixels.x)+1);
	Interpolate(leftPixels, rightPixels, rowPixels);
	for (int i=0; i<rowPixels.size(); i++)
	{
		PutPixelSDL(screen, rowPixels[i].x, rowPixels[i].y, currentColor);
	}
	/*
	for (int i=0; i<leftPixels.size(); i++)
	{
		int _x = leftPixels[i].x;
		int _y = leftPixels[i].y;

		for (int j=_x; j<=rightPixels[i].x; j++)
		{

			PutPixelSDL(screen, j, _y, currentColor);
		}
	}*/
}

void DrawPolygon(const vector<vec3>& vertices)
{
	int V = vertices.size();
	currentColor = vec3(vertices.color);

	vector<glm::ivec2> vertexPixels(V);
	for (int i=0; i<V; i++)
	{
		VertexShader(vertices[i], vertexPixels[i]);
	}

	vector<glm::ivec2> leftPixels;
	vector<glm::ivec2> rightPixels;
	ComputePolygonRows(vertexPixels, leftPixels, rightPixels);
	DrawPolygonRows(leftPixels, rightPixels);
}

void InterpolateP(Pixel a, Pixel b, vector<Pixel>& result)
{
	int N = result.size();
	vec3 step = vec3(b-a) / float(max(N-1,1));
	vec3 current(a);

	for (int i=0; i<N; ++i)
	{
		result[i] = current;
		current += step;
	}
}

void VertexShaderP(const vec3& v, Pixel p)
{
	// take 3D position of a vertex v and compute its 2D img position and store it in p
	vec3 pos = (v - cameraPos) * R;

	p.zinv = 1/pos.z;
	p.x = (focalLength * ( pos.x /  pos.z )) + (SCREEN_WIDTH / 2);
	p.y = (focalLength * ( pos.y /  pos.z )) + (SCREEN_HEIGHT / 2);
	
	depthBuffer[p.y][p.x] = p.zinv;
}

void ComputePolygonRowsP(const vector<Pixel>& vertexPixels,
			      vector<Pixel>& leftPixels,
			      vector<Pixel>& rightPixels)
{
	int minY = 999;
	int maxY = -999;
	int rows;

	// 1. Find max & min y-value of the polygon and compute the number of rows it occupies
	for (int i=0; i<vertexPixels.size(); i++)
	{
		minY = glm::min(minY, vertexPixels[i].y);
		maxY = glm::max(maxY, vertexPixels[i].x);
	}
	rows = maxY - minY + 1;

	// 2. Resize leftPixels and rightPixels so that they have an element for each row
	leftPixels.resize(rows);
	rightPixels.resize(rows);

	// 3. Initialize the x-coordinates in leftPixels to some really large value
	//    and the x-coordinates in the rightPixels to some really small values
	for (int i=0; i<rows; i++)
	{
		leftPixels[i].x = +numeric_limits<int>::max();
		rightPixels[i].x = -numeric_limits<int>::max();

		leftPixels[i].y = minY + i;
		rightPixels[i].y = minY + i;
	}

	// 4. Loop through all edges of the polygon and use linear interpolation to find
	//    the x-coordinate for each row it occupies. Update the corresponding values in rightPixels
	//    and leftPixels
	for (int i=0; i<vertexPixels.size(); i++)
	{
		int j = (i+1)%V; //next vertex
		Pixel delta = glm::abs(vertexPixels[i]-vertexPixels[j]);
		int pixels = glm::max(delta.x, delta.y) + 1;

		vector<Pixel> line(pixels);
		InterpolateP(vertexPixels[i], vertexPixels[j], line);

		// check for every point in line
		for (int a=0; a<line.size(); a++)
		{
			// check for every point in left & right pixels
			for (int b=0; b<rows; b++)
			{
				// if they are in the same y line
				//how about the zinv?
				if (leftPixels[b].y == line[a].y)
				{
					int leftX = glm::min(leftPixels[b].x, line[a].x);
					int rightX = glm::max(rightPixels[b].x, line[a].x);

					leftPixels[b].x = leftX;
					rightPixels[b].x = rightX;
				}
			}
		}
	}
}

void DrawPolygonRowsP(const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels)
{
	vector<Pixel> rowPixels((rightPixels.x - leftPixels.x)+1);
	Interpolate(leftPixels, rightPixels, rowPixels);
	for (int i=0; i<rowPixels.size(); i++)
	{
		PutPixelSDL(screen, rowPixels[i].x, rowPixels[i].y, currentColor);
	}
	/*
	for (int i=0; i<leftPixels.size(); i++)
	{
		int _x = leftPixels[i].x;
		int _y = leftPixels[i].y;

		for (int j=_x, j<=rightPixels[i].x; j++)
		{

			PutPixelSDL(screen, j, _y, currentColor);
		}
	}
	*/
}

void PixelShader(const Pixel& p)
{
	int x = p.x;
	int y = p.y;
	if (p.zinv > depthBuffer[y][x])
	{
		depthBuffer[y][x] = p.zinv;
		PutPixelSDL(screen, x, y, currentColor);
	}
}