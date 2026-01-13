/**
 * @file Core/Image.hpp
 * @brief Simple 2D image buffer storing RGB-like pixels in a flat array.
 *
 * The Image class represents a 2D raster with row-major storage (y * width + x).
 * Pixels are stored as Vector3 (typically linear RGB).
 *
 * Notes:
 * - Accessors do not perform bounds checking for performance.
 * - Use isValid() to verify dimensions and buffer size consistency.
 */


#pragma once

#include <vector>
#include <utility>

#include "Core/Types.hpp"
#include "Core/Assert.hpp"

namespace rayt {

	/**
	 * @class Image
	 * @brief A class representing a 2D image buffer.
	 *
	 * This class stores image data as a flat vector of Vector3 pixels (typically RGB).
	 * It provides methods to access pixels using 2D coordinates (x, y).
	 */
	class Image {
	public:
		/**
		 * @brief Default constructor. Creates an empty image.
		 */
		Image() = default;

		/**
		 * @brief Constructs an image with specific dimensions and pixel data.
		 *
		 * @param w The width of the image.
		 * @param h The height of the image.
		 * @param pixels A vector containing the pixel data. Its size must be w * h.
		 */
		Image(int w, int h, std::vector<Vector3> pixels)
			:m_width(w), m_height(h), m_pixels(std::move(pixels))
		{
			Assert(isValid() && "Image: pixel buffer size must be width * height.");
		}

		/**
		 * @brief Checks if the image data is valid.
		 *
		 * Verifies that the dimensions are positive and that the number of pixels
		 * matches the expected size (width * height).
		 *
		 * @return True if the image is valid, false otherwise.
		 */
		[[nodiscard]] bool isValid() const noexcept {
			if (m_width <= 0 || m_height <= 0) return false;

			const auto expected = static_cast<std::size_t>(m_width) * static_cast<std::size_t>(m_height);
			return m_pixels.size() == expected;
		}

		/**
		 * @brief Gets the width of the image.
		 * @return The width in pixels.
		 */
		[[nodiscard]] int width() const noexcept { return m_width; }

		/**
		 * @brief Gets the height of the image.
		 * @return The height in pixels.
		 */
		[[nodiscard]] int height() const noexcept { return m_height; }

		/**
		 * @brief Accesses the pixel at the specified coordinates (read-only).
		 *
		 * The coordinates are mapped to the 1D vector index using: y * width + x.
		 * No bounds checking is performed for performance reasons.
		 *
		 * @param x The x-coordinate (column).
		 * @param y The y-coordinate (row).
		 * @return A const reference to the pixel at (x, y).
		 */
		[[nodiscard]] const Vector3& at(int x, int y) const noexcept {
			Assert(x >= 0 && x < m_width && y >= 0 && y < m_height);
			return m_pixels[static_cast<std::size_t>(y) * static_cast<std::size_t>(m_width) + static_cast<std::size_t>(x)];
		}

		/**
		 * @brief Accesses the pixel at the specified coordinates (modifiable).
		 *
		 * The coordinates are mapped to the 1D vector index using: y * width + x.
		 * No bounds checking is performed for performance reasons.
		 *
		 * @param x The x-coordinate (column).
		 * @param y The y-coordinate (row).
		 * @return A reference to the pixel at (x, y).
		 */
		[[nodiscard]] Vector3& at(int x, int y) noexcept {
			return m_pixels[static_cast<std::size_t>(y) * static_cast<std::size_t>(m_width) + static_cast<std::size_t>(x)];
		}

		/**
		 * @brief Gets the raw pixel data (read-only).
		 * @return A const reference to the vector of pixels.
		 */
		[[nodiscard]] const std::vector<Vector3>& pixels() const noexcept { return m_pixels; }

		/**
		 * @brief Gets the raw pixel data (modifiable).
		 * @return A reference to the vector of pixels.
		 */
		[[nodiscard]] std::vector<Vector3>& pixels() noexcept { return m_pixels; }

		/// Pointer access (useful for IO).
		[[nodiscard]] const Vector3* data() const noexcept { return m_pixels.data(); }
		[[nodiscard]] Vector3* data() noexcept { return m_pixels.data(); }

	private:
		int m_width = 0;
		int m_height = 0;
		std::vector<Vector3> m_pixels;
	};
}