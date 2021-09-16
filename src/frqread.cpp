// MacRes v0.0.2   9/1/21
#include <algorithm>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "frqread.h"

#pragma warning(disable:4996)

/**
 * Helper function to convert a byte buffer to an integer. Supports both endians.
 */
int bufToInt(unsigned char *buf, char endian)
{
	int result;
	if (endian == 'l') 
	{
		result = buf[0] + (buf[1] << 8) + (buf[2] << 16) + (buf[3] << 24);
	}
	else
	{
		result = (buf[0] << 24) + (buf[1] << 16) + (buf[2] << 8) + buf[3];
	}
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
 * param sample_rate: sample rate of the original wav file.
 * param num_samples: number of samples read from the wav file.
 * param offset_ms: ms to start reading from.
 * param num_frames: put the number of frq values read here.
 * params ms_per_frq: put the ms per frq value (read from the file) here.
 * Returns the f0 curve read, or NULL on failure.
 */
double * ReadFrqFile(
	char *wav_filename,
	int sample_rate,
	int num_samples,
	int offset_ms,
	int *num_frames,
	double *ms_per_frq)
{
	char endian = 'l'; // 'l' for little-endian, 'b' for big-endian.
	int i; // For loops.
	int filename_len;
	int samples_per_frq;
	int total_num_frames, expected_num_frames, start_frame, end_frame;
	double ms_per_frame;
	double avg_frq;

	char *frq_filename;
	char header_buf[5]; // One more byte than needed.
	header_buf[4] = '\0'; // Needs to be null-terminated for strcmp to work.
	unsigned char int_buf[4];
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
	
	FILE *file;
	file = fopen(frq_filename, "rb");
	if (file == NULL)
	{
		printf("No .frq file found, generating f0.\n");
		return NULL;
	}
	printf("Reading f0 curve from .frq file.\n");
	
	// Check headers.
	size_t read_result = fread(header_buf, sizeof(char), 4, file); // "FREQ"
	assert(read_result == 4);
	if (strcmp(header_buf, "FREQ") != 0)
	{
		fclose(file);
		free(frq_filename);
		fprintf(stderr, "Error: Missing FREQ header in input file.\n");
		return NULL;
	}
	read_result = fread(header_buf, sizeof(char), 4, file); // "0003"
	assert(read_result == 4);
	if (strcmp(header_buf, "0003") != 0)
	{
		fclose(file);
		free(frq_filename);
		fprintf(stderr, "Error: Missing 0003 header in input file.\n");
		return NULL;
	}
	read_result = fread(int_buf, sizeof(char), 4, file); // Samples per frq.
	assert(read_result = 4);
	if (int_buf[0] == 0 && int_buf[1] == 1 && int_buf[2] == 0 && int_buf[3] == 0)
	{
		endian = 'l';
	}
	else if (int_buf[0] == 0 && int_buf[1] == 0 && int_buf[2] == 1 && int_buf[3] == 0)
	{
		endian = 'b';
	}
	else
	{
		fclose(file);
		free(frq_filename);
		fprintf(stderr, "Error: Only supports samples_per_frq value of 256.\n");
		return NULL;
	}
	samples_per_frq = bufToInt(int_buf, endian);
	ms_per_frame = samples_per_frq * 1000.0 / sample_rate;
	// Average frequency.
	read_result = fread(double_buf, sizeof(char), 8, file);
	assert(read_result = 8);
	avg_frq = bufToDouble(double_buf);
	// Empty space.
	fseek(file, 16, SEEK_CUR);
	// Number of frames.
	read_result = fread(int_buf, sizeof(char), 4, file);
	assert(read_result = 4);
	total_num_frames = bufToInt(int_buf, endian);
	// Total samples divided by samples_per_frq, rounded up.
	expected_num_frames = (int)((double)num_samples / (double)samples_per_frq) + 1;
	start_frame = (int) ((double)offset_ms / ms_per_frame);
	end_frame = start_frame + expected_num_frames;
	if (end_frame > total_num_frames)
	{
		fclose(file);
		free(frq_filename);
		fprintf(stderr, "Error: Need %d frames, found %d.\n", end_frame, total_num_frames);
		return NULL;
	}
	
	size_t total_num_bytes = sizeof(double) * total_num_frames * 2;
	char *f0_buf = (char *) malloc(sizeof(char) * total_num_bytes);
	double *f0 = (double *)malloc(sizeof(double) * expected_num_frames);
	
	// Read everything else from the file into memory.
	read_result = fread(f0_buf, sizeof(char), total_num_bytes, file);
	assert(read_result == total_num_bytes);
	
	// Extract every other double as the f0 curve.
	int cur_byte;
	for (int cur_frame = start_frame; cur_frame < end_frame; cur_frame++)
	{
		cur_byte = cur_frame * 16;
		f0[cur_frame] = bufToDouble(f0_buf + cur_byte);
	}
	
	// Clean up.
	free(frq_filename);
	free(f0_buf);
	fclose(file);
	
	// Populate return values.
	*num_frames = expected_num_frames;
	*ms_per_frq = ms_per_frame;
	return f0;
}
