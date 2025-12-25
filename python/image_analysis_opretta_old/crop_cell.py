

def crop_cell(
        mask: np.ndarray, # HW
        img: np.ndarray, # HWC, HW, output will preserve original ndim
        target_size: int = None,
        pad_square: bool = True,
        rotate: chr = None,
        exclude_border: bool = True, 
        max_cell_per_image: int = None,
        ) -> list[list[int], list[np.ndarray]]:
    """
    Crop cell
    Args:
        mask: np.ndarray, shape HW
        img: np.ndarray, shape HWC, HW, output will preserve original ndim
        target_size: int, target size
        pad_square: bool, pad square
        rotate: chr, rotate direction
        exclude_border: bool, exclude border
        max_cell_per_image: int, max cell per image
    Returns:
        list[list[int], list[np.ndarray]], cell ids and cropped cells
    """
    # check input
    if img.ndim == 2:
        raw_ndim = 2 # record original ndim
        img = np.expand_dims(img, 2) # HWC
    assert img.shape[:2] == mask.shape, "mask and img have different shape"
    
    # Clear border cell
    if exclude_border:
        mask = clear_border(mask)
    
    cell_ids = np.unique(mask)
    cell_ids = cell_ids[cell_ids != 0].ravel().tolist() # exclude background 0

    # Sample cell
    if max_cell_per_image and len(cell_ids) > max_cell_per_image:
        cell_ids = sample(cell_ids, max_cell_per_image)

    cropped_cells =[]
    for cell_id in cell_ids:
        cell_mask = (mask == cell_id).astype(np.uint8)
        
        # Get bounding box of the cell
        coords = np.argwhere(cell_mask)
        y_min, x_min = coords.min(axis=0)
        y_max, x_max = coords.max(axis=0)
        
        # Crop
        cropped_mask = cell_mask[y_min:y_max+1, x_min:x_max+1]
        cropped_cell = img[y_min:y_max+1, x_min:x_max+1, :] * cropped_mask[:,:,None]
        # imshow(cropped_cell[0,:,:])
        
        # Rotate mask and cell
        if rotate:
            cropped_mask, cropped_cell = rotate_cell(
                cropped_mask, cropped_cell, direction=rotate)
            # update coords
            coords = np.argwhere(cropped_mask)
            y_min, x_min = coords.min(axis=0)
            y_max, x_max = coords.max(axis=0)
        
        if pad_square:
            # Calculate size of the square
            height = y_max - y_min + 1
            width = x_max - x_min + 1
            square_size = max(height, width)
            # Calculate padding
            pad_y = (square_size - height) // 2
            pad_x = (square_size - width) // 2
            
            # Pad to square
            cropped_cell = np.pad(
                cropped_cell, ((pad_y, square_size - height - pad_y), (pad_x, square_size - width - pad_x), (0, 0)),
                mode='constant', constant_values=0)
            # cropped_mask = np.pad(
            #     cropped_mask, ((pad_y, square_size - height - pad_y), (pad_x, square_size - width - pad_x)),
            #     mode='constant', constant_values=0)
        
        # Resize to target size
        if target_size:
            cropped_cell = resize(
                cropped_cell, (target_size, target_size, cropped_cell.shape[2]), 
                preserve_range=True, anti_aliasing=True)
            # # Smooth resized mask edge
            # smoothed_mask = gaussian(cropped_padded_resized[0,:,:], sigma=1)
            # imshow(smoothed_mask)
            # smoothed_mask = (smoothed_mask > 0.5).astype(int)
            # imshow(smoothed_mask)
            # cropped_padded_resized[0,:,:] = smoothed_mask
            # imshow(cropped_padded_resized[0, :, :])
        
        # restore ndim
        if 'raw_ndim' in locals():
            cropped_cell = cropped_cell[:,:,0]

        cropped_cells.append(cropped_cell)

    return [cell_ids, cropped_cells]



def rotate_cell(
        cropped_mask: np.ndarray,
        cropped_cell: np.ndarray,
        direction: str = "horizonal"
        ) -> list[np.ndarray, np.ndarray]:
    """
    Rotate cell
    Args:
        cropped_mask: np.ndarray, shape HW
        cropped_cell: np.ndarray, shape HWC, HW
        direction: str, direction to rotate
    Returns:
        list[np.ndarray, np.ndarray], rotated mask and cell
    """
    # check input
    assert cropped_cell.shape[:2] == cropped_mask.shape, "mask and cell have different shape"

    # Get object properties
    props = regionprops(cropped_mask)
    
    # For one object, get its orientation
    orientation = props[0].orientation  # Orientation in radians

    # Calculate padding to prevent cropping during rotation
    max_dim = max(cropped_cell.shape)
    pad_width = max_dim
    if cropped_cell.ndim == 3:
        padded_img = np.pad(cropped_cell, ((pad_width, pad_width), (pad_width, pad_width), (0,0)), mode='constant', constant_values=0)
    else:
        padded_img = np.pad(cropped_cell, ((pad_width, pad_width), (pad_width, pad_width)), mode='constant', constant_values=0)
    padded_mask = np.pad(cropped_mask, ((pad_width, pad_width), (pad_width, pad_width)), mode='constant', constant_values=0)

    # Rotate image to align long axis horizontally
    angle_degrees = -np.degrees(orientation)
    if direction == "horizonal":
        angle_degrees = angle_degrees + 90
    rotated_img = ndi.rotate(padded_img, angle_degrees, reshape=False, mode='nearest', axes=(0, 1))
    rotated_mask = ndi.rotate(padded_mask, angle_degrees, reshape=False, mode='constant', cval=0)

    # Crop back to original image dimensions
    coords = np.argwhere(rotated_mask)
    y_min, x_min = coords.min(axis=0)
    y_max, x_max = coords.max(axis=0)

    if cropped_cell.ndim == 3:
        rotated_img = rotated_img[y_min:y_max+1, x_min:x_max+1, :]
    else:
        rotated_img = rotated_img[y_min:y_max+1, x_min:x_max+1]
    rotated_mask = rotated_mask[y_min:y_max+1, x_min:x_max+1]

    return [rotated_mask, rotated_img]

