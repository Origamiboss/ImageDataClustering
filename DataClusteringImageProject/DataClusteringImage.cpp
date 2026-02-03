#include <chrono>
#include <climits>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <fstream>
#include <math.h>
#include <string.h>
#include <cstdlib>
#include <time.h>
#include <float.h>

using namespace std;

typedef unsigned char uchar;
typedef unsigned long ulong;

typedef struct
{
	double red, green, blue;
} RGB_Pixel;

typedef struct
{
	int width, height;
	int size;
	RGB_Pixel* data;
} RGB_Image;
typedef struct {
	double* d1;
	double* d2;
	int* nearest_medoid;
} PAM_Workspace;
typedef struct
{
	int size;
	int SSE;
	double dist_to_nearest_medoid;
	double dist_to_second_nearest_medoid;
	int medoid_Index;
} RGB_Cluster;

/* Mersenne Twister related constants */
#define N 624
#define M 397
#define MAXBIT 30
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */
#define MAX_RGB_DIST 195075
#define NUM_RUNS 100

static ulong mt[N]; /* the array for the state vector  */
static int mti = N + 1; /* mti == N + 1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */
void init_genrand(ulong s)
{
	mt[0] = s & 0xffffffffUL;
	for (mti = 1; mti < N; mti++)
	{
		mt[mti] =
			(1812433253UL * (mt[mti - 1] ^ (mt[mti - 1] >> 30)) + mti);
		/* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
		/* In the previous versions, MSBs of the seed affect   */
		/* only MSBs of the array mt[].                        */
		/* 2002/01/09 modified by Makoto Matsumoto             */
		mt[mti] &= 0xffffffffUL;
		/* for >32 bit machines */
	}
}

ulong genrand_int32(void)
{
	ulong y;
	static ulong mag01[2] = { 0x0UL, MATRIX_A };
	/* mag01[x] = x * MATRIX_A  for x = 0, 1 */

	if (mti >= N)
	{ /* generate N words at one time */
		int kk;

		if (mti == N + 1)
		{
			/* if init_genrand ( ) has not been called, */
			init_genrand(5489UL); /* a default initial seed is used */
		}

		for (kk = 0; kk < N - M; kk++)
		{
			y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
			mt[kk] = mt[kk + M] ^ (y >> 1) ^ mag01[y & 0x1UL];
		}

		for (; kk < N - 1; kk++)
		{
			y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
			mt[kk] = mt[kk + (M - N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
		}

		y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
		mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ mag01[y & 0x1UL];
		mti = 0;
	}

	y = mt[mti++];

	/* Tempering */
	y ^= (y >> 11);
	y ^= (y << 7) & 0x9d2c5680UL;
	y ^= (y << 15) & 0xefc60000UL;
	y ^= (y >> 18);

	return y;
}

double genrand_real2(void)
{
	return genrand_int32() * (1.0 / 4294967296.0);
	/* divided by 2^32 */
}

/* Function for generating a bounded random integer between 0 and RANGE */
/* Source: http://www.pcg-random.org/posts/bounded-rands.html */

uint32_t bounded_rand(const uint32_t range)
{
	uint32_t x = genrand_int32();
	uint64_t m = ((uint64_t)x) * ((uint64_t)range);
	uint32_t l = (uint32_t)m;

	if (l < range)
	{
		//Without the (int32_t) The compilers flags this as an error because you are trying to make an unsigned int negative
		uint32_t t = -(int32_t)range;

		if (t >= range)
		{
			t -= range;
			if (t >= range)
			{
				t %= range;
			}
		}

		while (l < t)
		{
			x = genrand_int32();
			m = ((uint64_t)x) * ((uint64_t)range);
			l = (uint32_t)m;
		}
	}

	return m >> 32;
}

RGB_Image* read_PPM(const char* filename)
{
	uchar byte;
	char buff[16];
	int c, max_rgb_val, i = 0;
	FILE* fp;
	RGB_Pixel* pixel;
	RGB_Image* img;

	fopen_s(&fp, filename, "rb");
	if (!fp)
	{
		fprintf(stderr, "Unable to open file '%s'!\n", filename);
		exit(EXIT_FAILURE);
	}

	/* read image format */
	if (!fgets(buff, sizeof(buff), fp))
	{
		perror(filename);
		exit(EXIT_FAILURE);
	}

	/*check the image format to make sure that it is binary */
	if (buff[0] != 'P' || buff[1] != '6')
	{
		fprintf(stderr, "Invalid image format (must be 'P6')!\n");
		exit(EXIT_FAILURE);
	}

	img = (RGB_Image*)malloc(sizeof(RGB_Image));
	if (!img)
	{
		fprintf(stderr, "Unable to allocate memory!\n");
		exit(EXIT_FAILURE);
	}

	/* skip comments */
	c = getc(fp);
	while (c == '#')
	{
		while (getc(fp) != '\n');
		c = getc(fp);
	}

	ungetc(c, fp);

	/* read image dimensions */
	if (fscanf_s(fp, "%u %u", &img->width, &img->height) != 2)
	{
		fprintf(stderr, "Invalid image dimensions ('%s')!\n", filename);
		exit(EXIT_FAILURE);
	}

	/* read maximum component value */
	if (fscanf_s(fp, "%d", &max_rgb_val) != 1)
	{
		fprintf(stderr, "Invalid maximum R, G, B value ('%s')!\n", filename);
		exit(EXIT_FAILURE);
	}

	/* validate maximum component value */
	if (max_rgb_val != 255)
	{
		fprintf(stderr, "'%s' is not a 24-bit image!\n", filename);
		exit(EXIT_FAILURE);
	}

	while (fgetc(fp) != '\n');

	/* allocate memory for pixel data */
	img->size = img->height * img->width;
	img->data = (RGB_Pixel*)malloc(img->size * sizeof(RGB_Pixel));

	if (!img)
	{
		fprintf(stderr, "Unable to allocate memory!\n");
		exit(EXIT_FAILURE);
	}

	/* Read in pixels using buffer */
	while (fread(&byte, 1, 1, fp) && i < img->size)
	{
		pixel = &img->data[i];
		pixel->red = byte;
		fread(&byte, 1, 1, fp);
		pixel->green = byte;
		fread(&byte, 1, 1, fp);
		pixel->blue = byte;
		i++;
	}

	fclose(fp);

	return img;
}

void write_PPM(const RGB_Image* img, const char* filename)
{
	uchar byte;
	FILE* fp;

	fopen_s(&fp, filename, "wb");
	if (!fp)
	{
		fprintf(stderr, "Unable to open file '%s'!\n", filename);
		exit(EXIT_FAILURE);
	}

	fprintf(fp, "P6\n");
	fprintf(fp, "%d %d\n", img->width, img->height);
	fprintf(fp, "%d\n", 255);

	for (int i = 0; i < img->size; i++)
	{
		byte = (uchar)img->data[i].red;
		fwrite(&byte, sizeof(uchar), 1, fp);
		byte = (uchar)img->data[i].green;
		fwrite(&byte, sizeof(uchar), 1, fp);
		byte = (uchar)img->data[i].blue;
		fwrite(&byte, sizeof(uchar), 1, fp);
	}

	fclose(fp);
}

/* Function to generate random cluster centers. */
RGB_Cluster* gen_rand_centers(const RGB_Image* img, PAM_Workspace* ws, const int k) {
	RGB_Pixel rand_pixel;
	RGB_Cluster* cluster;
	cluster = (RGB_Cluster*)malloc(k * sizeof(RGB_Cluster));

	bool* used = (bool*)calloc(img->size, sizeof(bool));
	for (int i = 0; i < k; i++) {
		/* Make the initial guesses for the centers, m1, m2, ..., mk */
		//Make sure that the same pixel is not chosen twice
		int idx;
		do {
			idx = bounded_rand(img->size);
		} while (used[idx]);

		used[idx] = true;
		//set the index
		cluster[i].medoid_Index = idx;

		
		
		// cout << "\nCluster Centers: " << cluster[i].center.red << ", " << cluster[i].center.green <<", "<<  cluster[i].center.blue;
	}
	/* Set the number of points assigned to k cluster to zero, n1, n2, ..., nk */
	for (int i = 0; i < k; i++) {
		cluster[i].size = 0;
		cluster[i].SSE = 0;
		cluster[i].dist_to_nearest_medoid = DBL_MAX;
		cluster[i].dist_to_second_nearest_medoid = DBL_MAX;
	}
	//free up used
	free(used);
	//Set each pixel index
	ws->d1 = (double*)malloc(img->size * sizeof(double));
	ws->d2 = (double*)malloc(img->size * sizeof(double));
	ws->nearest_medoid = (int*)malloc(img->size * sizeof(int));
	for (int i = 0; i < img->size; i++) {
		ws->d1[i] = DBL_MAX;
		ws->d2[i] = DBL_MAX;

		for (int m = 0; m < k; m++) {
			int medoidIdx = cluster[m].medoid_Index;

			double dr = img->data[i].red - img->data[medoidIdx].red;
			double dg = img->data[i].green - img->data[medoidIdx].green;
			double db = img->data[i].blue - img->data[medoidIdx].blue;

			double dist = dr * dr + dg * dg + db * db;

			if (dist < ws->d1[i]) {
				ws->d2[i] = ws->d1[i];
				ws->d1[i] = dist;
				ws->nearest_medoid[i] = m;
			}
			else if (dist < ws->d2[i]) {
				ws->d2[i] = dist;
			}
		}
	}

	return(cluster);
}

int* get_rand_batch(const RGB_Image* img, const int batch_size) {
	int* batch = (int*)malloc(sizeof(int) * batch_size);
	if (!batch) return NULL;

	for (int i = 0; i < batch_size; i++) {
		batch[i] = bounded_rand(img->size);
	}

	return batch;
}


//Recolors the image to match the clusters
RGB_Image* map_pixels(const RGB_Image* src, RGB_Cluster* clusters, const int num_colors) {
	RGB_Image* dst = (RGB_Image*)malloc(sizeof(RGB_Image));
	//Make my new image
	dst->width = src->width;
	dst->height = src->height;
	dst->size = src->size;
	dst->data = (RGB_Pixel*)malloc(sizeof(RGB_Pixel) * dst->size);

	for (int i = 0; i < src->size; i++) {
		int closest = 0;

		// Find the closest cluster center
		// Initialize with the max distance
		double closestDist = DBL_MAX;

		for (int c = 1; c < num_colors; c++) {
			//Get the squared distance
			double red_diff = src->data[i].red - src->data[clusters[c].medoid_Index].red;
			double green_diff = src->data[i].green - src->data[clusters[c].medoid_Index].green;
			double blue_diff = src->data[i].blue - src->data[clusters[c].medoid_Index].blue;
			double dist = red_diff * red_diff + green_diff * green_diff + blue_diff * blue_diff;
			if (dist < closestDist) {
				closestDist = dist;
				closest = c;
			}
		}

		// Assign the pixel to the cluster center color
		dst->data[i] = src->data[clusters[closest].medoid_Index];
	}
	return dst;
}
/*
   For application of the batchk k-means algorithm to color quantization, see
   M. E. Celebi, Improving the Performance of K-Means for Color Quantization,
   Image and Vision Computing, vol. 29, no. 4, pp. 260-271, 2011.
 */
 /* Color quantization using the batch k-means algorithm */

// I have copied my k-means implementation here and have modified it to work with RGB images
void batch_kmedoids(const RGB_Image* img, PAM_Workspace* ws, const int num_colors, const int max_iters, RGB_Cluster* clusters)
{
	const int sizeOfInstance = 3; // RGB has 3 dimensions
	const int sizeOfBatch = 1000; // Batch size
	const double conversionThreshold = 0.001; // Convergence threshold

	//We will use Recolor Image later to assign the colors to the malloc image

	//i is the iteration we are on
	for (int i = 1; i <= max_iters; i++) {
		//Get a random batch that we will be working this for this iteration
		int* batch = (int*)malloc(sizeof(int) * sizeOfBatch);
		for (int i = 0; i < sizeOfBatch; i++) {
			batch[i] = bounded_rand(img->size);
		}
		//Decide on what we are swapping this iteration
		int changingCluster = bounded_rand(num_colors);
		int candidate = bounded_rand(img->size);
		double gain = 0.0;

		for (int b = 0; b < sizeOfBatch; b++) {
			int i = batch[b];

			int old_m = ws->nearest_medoid[i];

			double dr = img->data[i].red - img->data[candidate].red;
			double dg = img->data[i].green - img->data[candidate].green;
			double db = img->data[i].blue - img->data[candidate].blue;

			double dist_io = dr * dr + dg * dg + db * db;

			if (old_m == changingCluster) {
				// current medoid is being replaced
				gain += fmin(dist_io, ws->d2[i]) - ws->d1[i];
			}
			else {
				// unaffected medoid
				if (dist_io < ws->d1[i]) {
					gain += dist_io - ws->d1[i];
				}
			}
		}
		if (gain < 0.0) {
			clusters[changingCluster].medoid_Index = candidate;

			// update workspace
			for (int i = 0; i < img->size; i++) {
				double dr = img->data[i].red - img->data[candidate].red;
				double dg = img->data[i].green - img->data[candidate].green;
				double db = img->data[i].blue - img->data[candidate].blue;

				double dist_io = dr * dr + dg * dg + db * db;

				if (dist_io < ws->d1[i]) {
					ws->d2[i] = ws->d1[i];
					ws->d1[i] = dist_io;
					ws->nearest_medoid[i] = changingCluster;
				}
				else if (dist_io < ws->d2[i]) {
					ws->d2[i] = dist_io;
				}
			}
		}
		free(batch);
		cout << "Iteration: " << i << endl;
	}
	
}

void free_img(const RGB_Image* img) {
	/* Free Image Data*/
	free(img->data);

	/* Free Image Pointer*/
	delete(img);
}

int main(int argc, char* argv[])
{
	char* filename;						/* Filename Pointer*/
	int k;								/* Number of clusters*/
	//int batch_size;
	RGB_Image* img;
	RGB_Image* out_img;
	RGB_Cluster* cluster;

	if (argc == 3) {
		/* Image filename */
		filename = argv[1];

		/* k, number of clusters */
		k = atoi(argv[2]);

		//batch_size = atoi(argv[3]);
	}
	else if (argc > 3) {
		printf("Too many arguments supplied.\n");
		printf("DataClusteringImage.exe <file_name> <num_of_clusters>\n");
		//printf("DataClusteringImage.exe <file_name> <num_of_clusters> <batch_size>\n");
		return 0;
	}
	else {
		printf("Two arguments expected: image filename and number of clusters.\n");
		printf("DataClusteringImage.exe <file_name> <num_of_clusters>\n");
		//printf("DataClusteringImage.exe <file_name> <num_of_clusters> <batch_size>\n");
		printf("%d\n", argc);
		return 0;
	}

	srand(time(NULL));

	/* Print Args*/
	printf("%s %d\n", filename, k);

	/* Read Image*/
	img = read_PPM(filename);

	/* Test Batch K-Means*/
	/* Start Timer*/
	auto start = std::chrono::high_resolution_clock::now();

	/* Initialize centers */
	PAM_Workspace* importantPAMData = (PAM_Workspace*)malloc(sizeof(PAM_Workspace));
	cluster = gen_rand_centers(img, importantPAMData, k);

	/* Execute Batch K-means*/
	//RGB_Image* cluster_img = batch_kmeans(img, batch_size, k, INT_MAX, cluster);
	const int max_iters = 500;
	batch_kmedoids(img, importantPAMData, k, max_iters, cluster);

	//Now get the image based on the new clusters
	RGB_Image* cluster_img = map_pixels(img, cluster, k);
	
	
	

	/* Stop Timer*/
	auto stop = std::chrono::high_resolution_clock::now();

	/* Execution Time*/
	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

	free(cluster);
	printf("Execution Time (ms): %lld\n", elapsed.count());

	//Write the new file
	const char* outputFilename = "output.ppm";
	write_PPM(cluster_img, outputFilename);

	free_img(img);
	free_img(cluster_img);
	return 0;
}