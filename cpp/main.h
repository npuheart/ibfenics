/// @date 2025-01-01
/// @file main.h
/// @author Ma Pengfei (code@pengfeima.cn)
/// @version 0.1
/// @copyright Copyright (c) 2025 Ma Pengfei
///
/// @brief
///
///

#pragma once

// 排列方式：
// p1_x, p1_y, p2_x, p2_y, p3_x, p3_y

template <size_t DIM>
void linearInterpolation(double *points_local, const double *points_ref, const double *vertices, size_t num_points)
{
    for (std::size_t i = 0; i < num_points; ++i)
    {
        double alpha = 1.0 - points_ref[DIM * i] - points_ref[DIM * i + 1];
        double beta = points_ref[DIM * i];
        double gamma = points_ref[DIM * i + 1];

        for (size_t j = 0; j < DIM; j++)
        {

            points_local[DIM * i + j] = alpha * vertices[j] + beta * vertices[DIM + j] + gamma * vertices[2 * DIM + j];
            // printf("alpha: %f, beta: %f, gamma: %f\n", alpha, beta, gamma);
        }
    }
}
