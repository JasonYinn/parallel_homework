import numpy as np
import os

def get_rays_np(H, W, K, c2w):
    i, j = np.meshgrid(np.arange(W, dtype=np.float32), np.arange(H, dtype=np.float32), indexing='xy')
    dirs = np.stack([(i-K[0][2])/K[0][0], -(j-K[1][2])/K[1][1], -np.ones_like(i)], -1)
    # Rotate ray directions from camera frame to the world frame
    print ((dirs[..., np.newaxis, :] * c2w[:3,:3]).shape)
    p = dirs[0, 0]
    tmp = c2w 
    test = np.zeros((4,))
    test[:3] = p
    a = [0., 0., 0., 0.]
    for t in range(4):
        print (test[t])
        print (test[t] * c2w[:, t])
        # a = a + test[t] * c2w[:, t]
    # print (a)
    test = np.stack([test * tmp[0], test * tmp[1], test * tmp[2], test * tmp[3]], axis=-1)
    
    # for t in range(4):
        # print (c2w[:, t])
        # print (test[:, t])
    # for t in test:
    #     print (t)
    # print (test)
    test = test.sum(0)
    # print (test)
    rays_d = np.sum(dirs[..., np.newaxis, :] * c2w[:3,:3], -1)  # dot product, equals to: [c2w.dot(dir) for dir in dirs]
    # Translate camera frame's origin to the world frame. It is the origin of all rays.
    rays_o = np.broadcast_to(c2w[:3,-1], np.shape(rays_d))
    return rays_o, rays_d, test

W, H = 800, 800
camera_angle_x = 0.9807014465332031
focal = (.5 * W / np.tan(.5 * camera_angle_x))

W, H = W // 2, H // 2
focal = focal / 2.0
c2w = np.array([
    [0.9923269748687744, 0.04875849187374115, -0.11362089961767197, -0.3664160668849945], 
    [-0.1236410066485405, 0.39132946729660034, -0.911906898021698, -2.359870433807373], 
    [0.0, 0.9189580678939819, 0.39435532689094543, 1.2530806064605713], 
    [0.0, 0.0, 0.0, 1.0]])

K = np.array([
    [focal, 0, 0.5*W],
    [0, focal, 0.5*H],
    [0, 0, 1]
])
rays_o, rays_d, test = get_rays_np(H, W, K, c2w)
print (rays_d[0, 0])