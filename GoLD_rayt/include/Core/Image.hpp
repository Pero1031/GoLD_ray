#pragma once

#include <vector>

#include "Core/Types.hpp"
#include "Core/Constants.hpp"

namespace rayt{

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
			:m_width(w), m_height(h), m_pixels(std::move(pixels)){}

		/**
		 * @brief Checks if the image data is valid.
		 *
		 * Verifies that the dimensions are positive and that the number of pixels
		 * matches the expected size (width * height).
		 *
		 * @return True if the image is valid, false otherwise.
		 */
		bool isValid() const {
			return m_width > 0 && m_height > 0 &&
				static_cast<int>(m_pixels.size()) == m_width * m_height;
		}

		/**
		 * @brief Gets the width of the image.
		 * @return The width in pixels.
		 */
		int width() const { return m_width; }

		/**
		 * @brief Gets the height of the image.
		 * @return The height in pixels.
		 */
		int height() const { return m_height; }

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
		const Vector3& at(int x, int y) const {
			return m_pixels[y * m_width + x];
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
		Vector3& at(int x, int y) {
			return m_pixels[y * m_width + x];
		}

		/**
		 * @brief Gets the raw pixel data (read-only).
		 * @return A const reference to the vector of pixels.
		 */
		const std::vector<Vector3>& pixels() const { return m_pixels; }

		/**
		 * @brief Gets the raw pixel data (modifiable).
		 * @return A reference to the vector of pixels.
		 */
		std::vector<Vector3>& pixels() { return m_pixels; }

	private:
		int m_width = 0;
		int m_height = 0;
		std::vector<Vector3> m_pixels;
	};
}