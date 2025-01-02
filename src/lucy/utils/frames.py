''' frames.py  The original purpose of this program was to extract and randomly send frames of animal cruelty footage from cd ../..
    Copyright (C) 2024  github.com/brandongrahamcobb

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
'''
from utils.setup_logging import logger

import cv2
import os
import random

def extract_random_frames(video_path, output_dir, num_frames=1):
    try:
        logger.info(f'Starting to extract random frames from video: {video_path}')

        # Ensure the output directory exists
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
            logger.info(f'Created output directory: {output_dir}')

        # Open the video file
        cap = cv2.VideoCapture(video_path)
        if not cap.isOpened():
            logger.error(f'Failed to open video file: {video_path}')
            return []

        frame_count = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
        logger.debug(f'Total number of frames in video: {frame_count}')

        # Generate random frame indices
        frame_indices = random.sample(range(frame_count), num_frames)
        logger.debug(f'Selected random frame indices: {frame_indices}')

        extracted_frames = []

        # Extract and save frames
        for i, frame_idx in enumerate(frame_indices):
            logger.debug(f'Processing frame index: {frame_idx}')
            cap.set(cv2.CAP_PROP_POS_FRAMES, frame_idx)
            ret, frame = cap.read()

            if ret:
                frame_path = os.path.join(output_dir, f'frame_{i + 1}.jpg')
                cv2.imwrite(frame_path, frame)
                extracted_frames.append(frame_path)
                logger.info(f'Frame saved: {frame_path}')
            else:
                logger.warning(f'Failed to read frame at index: {frame_idx}')

        cap.release()
        logger.info(f'Completed frame extraction. Total frames extracted: {len(extracted_frames)}')

        return extracted_frames

    except Exception as e:
        logger.error(f'An error occurred while extracting frames: {e}')
        raise
