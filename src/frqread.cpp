// MacRes v0.0.2   9/1/21
#include <algorithm>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "frqread.h"

#pragma warning(disable:4996)

/**
 * Helper function to convert a byte buffer to an integer.
 */
int bufToInt(unsigned char *buf)
{
	int result = static_cast<int>(buf[0] << 24 | buf[1] << 16 | buf[2] << 8 | buf[3]);
	return result;
}

/**
 * Helper function to convert a byte buffer to a double.
 */
double bufToDouble(char *buf)
{
	double result;
	std::copy(buf, buf + sizeof(double), reinterpret_cast<char *>(&result));
	return result;
}

/**
 * Reads a signal from a FRQ0003 file.
 * param wav_filename: file name of wav file.
 * param num_frames: number of frq frames this file should have
 * f0: Put the f0 curve from the file here. Should already be malloc'ed.
 * Returns 0 on success, -1 on failure.
 */
int ReadFrqFile(
	char *wav_filename,
	int num_frames,
	double *f0)
{
	int i; // For loops.
	int filename_len;
	double avg_frq;

	char *frq_filename;
	char header_buf[5]; // One more byte than needed.
	header_buf[4] = '\0'; // Needs to be null-terminated for strcmp to work.
	unsigned char integer_buf[4];
	char double_buf[8]; // Assume double size is 8.
	
	// Convert wav filename into the equivalent frq filename.
	// Example: "a.wav" -> "a_wav.frq"
	filename_len = strlen(wav_filename);
	frq_filename = (char *) malloc(sizeof(char) * (filename_len + 5));
	for (i = 0; i < filename_len; i++)
	{
		frq_filename[i] = wav_filename[i];
	}
	frq_filename[filename_len - 4] = '_';
	frq_filename[filename_len] = '.';
	frq_filename[filename_len + 1] = 'f';
	frq_filename[filename_len + 2] = 'r';
	frq_filename[filename_len + 3] = 'q';
	frq_filename[filename_len + 4] = '\0';
	printf("%s\n", frq_filename);
	
	FILE *file;
	file = fopen(frq_filename, "rb");
	if (file == NULL)
	{
		printf("No .frq file found, generating f0.\n");
		return -1;
	}
	printf("Reading f0 curve from .frq file.\n");
	
	// Check headers.
	size_t read_result = fread(header_buf, sizeof(char), 4, file); // "FREQ"
	assert(read_result == 4);
	if (strcmp(header_buf, "FREQ") != 0)
	{
		fclose(file);
		fprintf(stderr, "Error: Missing FREQ header in input file.\n");
		return -1;
	}
	read_result = fread(header_buf, sizeof(char), 4, file); // "0003"
	assert(read_result == 4);
	if (strcmp(header_buf, "0003") != 0)
	{
		fclose(file);
		fprintf(stderr, "Error: Missing 0003 header in input file.\n");
		return -1;
	}
	read_result = fread(integer_buf, sizeof(char), 4, file); // Samples per frq.
	assert(read_result = 4);
	if (bufToInt(integer_buf) != 256)
	{
		fclose(file);
		fprintf(stderr, "Error: Only supports samples_per_frq value of 256.\n");
		return -1;
	}
	read_result = fread(double_buf, sizeof(char), 8, file); // Average frequency.
	assert(read_result = 8);
	avg_frq = bufToDouble(double_buf);
	fseek(file, 16, SEEK_CUR); // Empty space.
	read_result = fread(integer_buf, sizeof(char), 4, file); // Number of frames.
	assert(read_result = 4);
	if (bufToInt(integer_buf) != num_frames)
	{
		fclose(file);
		fprintf(
			stderr,
			"Error: num frames read from .frq file is %d, should be %d.\n",
			bufToInt(integer_buf),
			num_frames);
		return -1;
	}
	
	size_t total_num_bytes = sizeof(double) * num_frames * 2;
	char *f0_buf = (char *) malloc(sizeof(char) * total_num_bytes);
	
	// Read everything else from the file into memory.
	read_result = fread(f0_buf, sizeof(char), total_num_bytes, file);
	assert(read_result == total_num_bytes);
	
	// Extract every other double as the f0 curve.
	int cur_byte;
	for (int cur_frame = 0; cur_frame < num_frames; cur_frame++)
	{
		cur_byte = cur_frame * 16;
		f0[cur_frame] = bufToDouble(f0_buf + cur_byte);
	}
	
	// Clean up.
	free(frq_filename);
	free(f0_buf);
	fclose(file);
	return 0;
}
