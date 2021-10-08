#include "Animation.hpp"
#include <igl/opengl/glfw/Viewer.h>

Eigen::MatrixXi Animation::get_edge_matrix()
{
  int total_edges = 0;
  for (int i = 0; i < gm.cells.size(); i++)
  {
    total_edges += 4;
    if (gm.cells[i].type == RIGID)
    {
      total_edges += 2;
    }
    if (gm.cells[i].type == ACTIVE)
    {
      total_edges += 1;
    }
  }
  Eigen::MatrixXi E = Eigen::MatrixXi::Zero(total_edges, 2);

  int index = 0;
  for (int i = 0; i < gm.cells.size(); i++)
  {
    // Add four boundary edges
    E(index, 0) = gm.cells[i].vertices(0);
    E(index, 1) = gm.cells[i].vertices(1);
    index++;
    E(index, 0) = gm.cells[i].vertices(1);
    E(index, 1) = gm.cells[i].vertices(2);
    index++;
    E(index, 0) = gm.cells[i].vertices(2);
    E(index, 1) = gm.cells[i].vertices(3);
    index++;
    E(index, 0) = gm.cells[i].vertices(3);
    E(index, 1) = gm.cells[i].vertices(0);
    index++;
    if (gm.cells[i].type == RIGID)
    { // Diagonal edges
      E(index, 0) = gm.cells[i].vertices(1);
      E(index, 1) = gm.cells[i].vertices(3);
      index++;
      E(index, 0) = gm.cells[i].vertices(0);
      E(index, 1) = gm.cells[i].vertices(2);
      index++;
    }
    if (gm.cells[i].type == ACTIVE)
    { // Single diagonal edges
      E(index, 0) = gm.cells[i].vertices(1);
      E(index, 1) = gm.cells[i].vertices(3);
      index++;
    }
  }

  return E;
}

std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::MatrixXd, Eigen::MatrixXd> Animation::animate_path()
{
  frame = frame  % res.size(); // properly bound frame

  // Initialize Points
  int points_size = res[frame].points.size() + to_trace.size() * res.size();
  for (auto path : target_paths) {
    points_size += path.size();
  }
  Eigen::MatrixXd P = Eigen::MatrixXd::Zero(points_size, 3);

  // Initialize Edges
  int edges_size = grid_edges.rows() + to_trace.size() * (res.size() - 1);
  for (auto path : target_paths) {
    edges_size += path.size() - 1;
  }
  Eigen::MatrixXi E = Eigen::MatrixXi::Zero(edges_size, 2);

  // Points on grid
  for (int i = 0; i < res[frame].points.size(); i++)
  {
    P(i, 0) = res[frame].points[i](0);
    P(i, 1) = res[frame].points[i](1);
    P(i, 2) = 0;
  }

  // Edges on grid
  E.topRows(grid_edges.rows()) = grid_edges;

  // Traced paths
  for (int i = 0; i < to_trace.size(); i++)
  { // ith traced path
    for (int j = 0; j < res.size(); j++)
    { // at time j
      P(res[frame].points.size() + i * res.size() + j, 0) = res[j].points[to_trace[i]](0);
      P(res[frame].points.size() + i * res.size() + j, 1) = res[j].points[to_trace[i]](1);
      P(res[frame].points.size() + i * res.size() + j, 2) = 0;
      if (j > 0) {
        E(grid_edges.rows() + i * res.size() + (j-1), 0) = res[frame].points.size() + i * res.size() + j - 1;
        E(grid_edges.rows() + i * res.size() + (j-1), 1) = res[frame].points.size() + i * res.size() + j;
      }
    }
  }

  // Target paths
  int i = 0;
  int j = 0;
  for (auto path : target_paths)
  {
    bool first = true;
    for (auto point : path)
    {
      P(res[frame].points.size() + to_trace.size() * res.size() + i, 0) = point[0];
      P(res[frame].points.size() + to_trace.size() * res.size() + i, 1) = point[1];
      P(res[frame].points.size() + to_trace.size() * res.size() + i, 2) = 0;
      if (!first) {
        E(grid_edges.rows() + to_trace.size() * (res.size() - 1) + j, 0) = res[frame].points.size() + to_trace.size() * res.size() + i - 1;
        E(grid_edges.rows() + to_trace.size() * (res.size() - 1) + j, 1) = res[frame].points.size() + to_trace.size() * res.size() + i;
        j++;
      } else {
        first = false;
      }
      i++;
    }
  }

  // // Edges on traced paths
  // for (int i = 0; i < to_trace.size(); i++)
  // { // ith traced path
  //   for (int j = 0; j < res.size() - 1; j++)
  //   { // edges between time j and j+1
  //     E(grid_edges.rows() + i * res.size() + j, 0) = res[frame].points.size() + i * res.size() + j;
  //     E(grid_edges.rows() + i * res.size() + j, 1) = res[frame].points.size() + i * res.size() + j + 1;
  //   }
  // }

  // Point Colors
  Eigen::MatrixXd PC = Eigen::MatrixXd::Zero(res[frame].points.size(), 3);
  for (int i : gm.anchors)
  {
    PC.row(i) = Eigen::RowVector3d(1, 0, 0);
  }

  // Edge Colors
  Eigen::MatrixXd EC = Eigen::MatrixXd::Zero(edges_size, 3);
  EC.bottomLeftCorner(edges_size - grid_edges.rows(), 1).setOnes();
  EC.bottomLeftCorner(edges_size - grid_edges.rows() - to_trace.size() * (res.size() - 1), 1).setZero();
  EC.bottomRightCorner(edges_size - grid_edges.rows() - to_trace.size() * (res.size() - 1), 1).setOnes();

  frame = (frame + 1) % res.size(); // Increment frame

  return std::make_tuple(P, E, PC, EC);
}

void Animation::animate()
{
  // Initialize viewer
  igl::opengl::glfw::Viewer viewer;
  viewer.data().point_size = 20;
  viewer.data().set_face_based(true);
  viewer.core().orthographic = true;
  viewer.core().toggle(viewer.data().show_lines);
  viewer.core().is_animating = true;
  viewer.core().background_color.setOnes();

  int curr_in_frame = 0;

  // Animation Callback
  viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer &) -> bool
  {
    curr_in_frame++;
    if (curr_in_frame == rate)
    { // slow controls frame rate
      // Get points and edges to draw
      auto tup = animate_path();
      auto points = std::get<0>(tup);
      auto edges = std::get<1>(tup);
      auto pointColors = std::get<2>(tup);
      auto edgeColors = std::get<3>(tup);
      // update points and edges
      viewer.data().set_points(points.topRows(gm.points.size()), pointColors);
      viewer.data().set_edges(points, edges, edgeColors);
      curr_in_frame = 0;
    }

    return false;
  };

  viewer.launch();
}