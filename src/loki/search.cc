#include "loki/search.h"
#include <valhalla/midgard/linesegment2.h>
#include <valhalla/midgard/util.h>

#include <unordered_set>
#include <list>
#include <math.h>

using namespace valhalla::midgard;
using namespace valhalla::baldr;
using namespace valhalla::sif;
using namespace valhalla::loki;

namespace {
//the cutoff at which we will assume the input is too far away from civilisation to be
//worth correlating to the nearest graph elements
constexpr float SEARCH_CUTOFF = 35000.f;
//the number of edges to have visited within a region to identify such a reason as not an island
constexpr size_t CONNECTED_REGION = 100;
//maximum number of candidates to keep track of when evaluating edges for correlation
constexpr size_t MAX_CANDIDATES = 20;
//during edge correlation, if you end up < 5 meters from the beginning or end of the
//edge we just assume you were at that node and not actually along the edge
//we keep it small because point and click interfaces are more accurate than gps input
constexpr float NODE_SNAP = 5.f;
//during side of street computations we figured you're on the street if you are less than
//5 meters (16) feet from the centerline. this is actually pretty large (with accurate shape
//data for the roads it might want half that) but its better to assume on street than not
constexpr float SIDE_OF_STREET_SNAP = 5.f;
//if you are this far away from the edge we are considering and you set a heading we will
//ignore it because its not really useful at this distance from the geometry
constexpr float NO_HEADING = 30.f;
//how much of the shape should be sampled to get heading
constexpr float HEADING_SAMPLE = 30.f;
//cone width to use for cosine similarity comparisons for favoring heading
constexpr float ANGLE_WIDTH = 88.f;

//TODO: move this to midgard and test the crap out of it
//we are essentially estimating the angle of the tangent
//at a point along a discretised curve. we attempt to mostly
//use the shape coming into the point on the curve but if there
//isnt enough there we will use the shape coming out of the it
float tangent_angle(size_t index, const PointLL& point, const std::vector<PointLL>& shape, bool forward) {
  //depending on if we are going forward or backward we choose a different increment
  auto increment = forward ? -1 : 1;
  auto first_end = forward ? shape.cbegin() : shape.cend() - 1 ;
  auto second_end = forward ? shape.cend() - 1 : shape.cbegin();

  //u and v will be points we move along the shape until we have enough distance between them or run out of points

  //move backwards until we have enough or run out
  float remaining = HEADING_SAMPLE;
  auto u = point;
  auto i = shape.cbegin() + index + forward;
  while(remaining > 0 && i != first_end) {
    //move along and see how much distance that added
    i += increment;
    auto d = u.Distance(*i);
    //are we done yet?
    if(remaining <= d) {
      auto coef = remaining / d;
      u = u.AffineCombination(1 - coef, coef, *i);
      return u.Heading(point);
    }
    //next one
    u = *i;
    remaining -= d;
  }

  //move forwards until we have enough or run out
  auto v = point;
  i = shape.cbegin() + index + !forward;
  while(remaining > 0 && i != second_end) {
    //move along and see how much distance that added
    i -= increment;
    auto d = v.Distance(*i);
    //are we done yet?
    if(remaining <= d) {
      auto coef = remaining / d;
      v = v.AffineCombination(1 - coef, coef, *i);
      return u.Heading(v);
    }
    //next one
    v = *i;
    remaining -= d;
  }

  return u.Heading(v);
}

bool heading_filter(const DirectedEdge* edge, const std::unique_ptr<const EdgeInfo>& info,
  const std::tuple<PointLL, float, int>& point, boost::optional<int> heading) {
  //if its far enough away from the edge, the heading is pretty useless
  if(!heading || std::get<1>(point) > NO_HEADING)
    return false;

  //get the angle of the shape from this point
  auto angle = tangent_angle(std::get<2>(point), std::get<0>(point), info->shape(), edge->forward());
  //we want the closest distance between two angles which can be had
  //across 0 or between the two so we just need to know which is bigger
  if(*heading > angle)
    return std::min(*heading - angle, (360.f - *heading) + angle) > ANGLE_WIDTH;
  return std::min(angle - *heading, (360.f - angle) + *heading) > ANGLE_WIDTH;
}

PathLocation::SideOfStreet flip_side(const PathLocation::SideOfStreet side) {
  if(side != PathLocation::SideOfStreet::NONE)
    return side == PathLocation::SideOfStreet::LEFT ? PathLocation::SideOfStreet::RIGHT : PathLocation::SideOfStreet::LEFT;
  return side;
}

PathLocation::SideOfStreet get_side(const DirectedEdge* edge, const std::unique_ptr<const EdgeInfo>& info,
  const std::tuple<PointLL, float, int>& point, const PointLL& original){

  //its so close to the edge that its basically on the edge
  if(std::get<1>(point) < SIDE_OF_STREET_SNAP)
    return PathLocation::SideOfStreet::NONE;

  //if the projected point is way too close to the begin or end of the shape
  //TODO: if the original point is really far away side of street may also not make much sense..
  if(std::get<0>(point).Distance(info->shape().front()) < SIDE_OF_STREET_SNAP ||
     std::get<0>(point).Distance(info->shape().back()) < SIDE_OF_STREET_SNAP)
    return PathLocation::SideOfStreet::NONE;

  //get the side TODO: this can technically fail for longer segments..
  //to fix it we simply compute the plane formed by the triangle
  //through the center of the earth and the two shape points and test
  //whether the original point is above or below the plane (depending on winding)
  auto index = std::get<2>(point);
  LineSegment2<PointLL> segment(info->shape()[index], info->shape()[index + 1]);
  return (segment.IsLeft(original) > 0) == edge->forward()  ? PathLocation::SideOfStreet::LEFT : PathLocation::SideOfStreet::RIGHT;
}

const NodeInfo* get_end_node(GraphReader& reader, const DirectedEdge* edge) {
  //the node could be in another tile so we grab that
  const auto tile = reader.GetGraphTile(edge->endnode());
  return tile->node(edge->endnode());
}

GraphId get_opposing(const GraphId& edge_id, const GraphTile*& tile, GraphReader& reader) {
  const auto* edge = tile->directededge(edge_id);
  if(edge->leaves_tile())
    tile = reader.GetGraphTile(edge->endnode());
  if(tile == nullptr)
    return {};
  auto id = edge->endnode();
  auto* node = tile->node(edge->endnode());
  id.fields.id = tile->node(id)->edge_index() + edge->opp_index();
  return id;
}

struct candidate_t{
  const GraphTile* tile = nullptr;
  const DirectedEdge* edge = nullptr;
  GraphId id;
  std::unique_ptr<const EdgeInfo> edge_info;
  std::tuple<PointLL, float, int> point{{}, std::numeric_limits<float>::max(), 0};

  //sorting
  bool operator<(const candidate_t& other) const { return std::get<1>(point) > std::get<1>(other.point); }
  //defaulting
  candidate_t() {}
  //moveing
  candidate_t(candidate_t&& c):tile(std::move(c.tile)),edge(std::move(c.edge)),id(std::move(c.id)),edge_info(std::move(c.edge_info)),point(std::move(c.point)) {}
  //copying
  candidate_t(const candidate_t& c):tile(c.tile),edge(c.edge),id(c.id),point(c.point) {}
  //opposing
  candidate_t(const candidate_t& c, GraphReader& reader):tile(c.tile),id(c.id) { id = get_opposing(id, tile, reader); if(id.Is_Valid()) edge = tile->directededge(id); }
};

std::vector<PathLocation::PathEdge> correlate_node(GraphReader& reader, const Location& location, const EdgeFilter& edge_filter, const candidate_t& closest){
  std::vector<PathLocation::PathEdge> correlated;
  std::list<PathLocation::PathEdge> heading_filtered;
  //now that we have a node we can pass back all the edges leaving and entering it
  const auto* tile = closest.tile->id().tileid() != closest.id.tileid() ? reader.GetGraphTile(closest.edge->endnode()) : closest.tile;
  const auto* node = tile->node(closest.edge->endnode());
  const auto* start_edge = tile->directededge(node->edge_index());
  const auto* end_edge = start_edge + node->edge_count();
  auto closest_point = closest.point;
  //for each edge at this node
  for(const auto* edge = start_edge; edge < end_edge; ++edge) {
    //get some info about this edge and the opposing
    GraphId id = tile->id();
    id.fields.id = node->edge_index() + (edge - start_edge);
    const GraphTile* other_tile = tile;
    auto other_id = get_opposing(id, other_tile, reader);
    if(!other_id.Is_Valid())
      continue;
    const auto* other_edge = other_tile->directededge(other_id);
    auto info = tile->edgeinfo(edge->edgeinfo_offset());

    //do we want this edge
    if(edge_filter(edge) != 0.0f) {
      PathLocation::PathEdge path_edge{std::move(id), 0.f, node->latlng(), std::get<1>(closest_point), PathLocation::NONE};
      std::get<2>(closest_point) = edge->forward() ? 0 : info->shape().size() - 2;
      if(!heading_filter(edge, info, closest_point, location.heading_))
        correlated.push_back(std::move(path_edge));
      else
        heading_filtered.emplace_back(std::move(path_edge));
    }

    //do we want the evil twin
    if(edge_filter(other_edge) != 0.0f) {
      PathLocation::PathEdge path_edge{std::move(other_id), 1.f, node->latlng(), std::get<1>(closest_point), PathLocation::NONE};
      std::get<2>(closest_point) = other_edge->forward() ? 0 : info->shape().size() - 2;
      if(!heading_filter(other_edge, tile->edgeinfo(edge->edgeinfo_offset()), closest_point, location.heading_))
        correlated.push_back(std::move(path_edge));
      else
        heading_filtered.emplace_back(std::move(path_edge));
    }
  }

  //TODO: now that we return multiple results we need to score
  //these lower or maybe not return them
  //if we have nothing because of heading we'll just ignore it
  if(correlated.size() == 0 && heading_filtered.size())
    for(auto& path_edge : heading_filtered)
      correlated.push_back(std::move(path_edge));

  //give it back
  return correlated;
}

std::vector<PathLocation::PathEdge> correlate_edge(GraphReader& reader, const Location& location, const EdgeFilter& edge_filter, const candidate_t& closest) {
  //now that we have an edge we can pass back all the info about it
  std::vector<PathLocation::PathEdge> correlated;
  if(closest.edge != nullptr){
    //we need the ratio in the direction of the edge we are correlated to
    double partial_length = 0;
    for(size_t i = 0; i < std::get<2>(closest.point); ++i)
      partial_length += closest.edge_info->shape()[i].Distance(closest.edge_info->shape()[i + 1]);
    partial_length += closest.edge_info->shape()[std::get<2>(closest.point)].Distance(std::get<0>(closest.point));
    partial_length = std::min(partial_length, static_cast<double>(closest.edge->length()));
    float length_ratio = static_cast<float>(partial_length / static_cast<double>(closest.edge->length()));
    if(!closest.edge->forward())
      length_ratio = 1.f - length_ratio;
    //side of street
    auto side = get_side(closest.edge, closest.edge_info, closest.point, location.latlng_);
    //correlate the edge we found
    std::list<PathLocation::PathEdge> heading_filtered;
    if(heading_filter(closest.edge, closest.edge_info, closest.point, location.heading_))
      heading_filtered.emplace_back(closest.id, length_ratio, std::get<0>(closest.point), side);
    else
      correlated.push_back(PathLocation::PathEdge{closest.id, length_ratio, std::get<0>(closest.point), std::get<1>(closest.point), side});
    //correlate its evil twin
    const GraphTile* other_tile = closest.tile;
    auto opposing_edge_id = get_opposing(closest.id, other_tile, reader);
    if(opposing_edge_id.Is_Valid()) {
      const auto* other_edge = other_tile->directededge(opposing_edge_id);
      if(edge_filter(other_edge) != 0.0f) {
        if(heading_filter(other_edge, closest.edge_info, closest.point, location.heading_))
          heading_filtered.emplace_back(opposing_edge_id, 1 - length_ratio, std::get<0>(closest.point), flip_side(side));
        else
          correlated.push_back(PathLocation::PathEdge{opposing_edge_id, 1 - length_ratio, std::get<0>(closest.point), std::get<1>(closest.point), flip_side(side)});
    const GraphTile* other_tile;
    auto opposing_edge_id = reader.GetOpposingEdgeId(closest_edge_id, other_tile);
    // Protect against empty tile (e.g. for extracts from tile sets)
    if (other_tile != nullptr) {
      const DirectedEdge* other_edge;
      if(opposing_edge_id.Is_Valid() && (other_edge = other_tile->directededge(opposing_edge_id)) && edge_filter(other_edge) != 0.0f) {
        if(heading_filter(other_edge, closest_edge_info, closest_point, location.heading_))
          heading_filtered.emplace_back(opposing_edge_id, 1 - length_ratio, std::get<0>(closest_point), std::get<1>(closest_point), flip_side(side));
        else
          correlated.edges.push_back(PathLocation::PathEdge{opposing_edge_id, 1 - length_ratio, std::get<0>(closest_point), std::get<1>(closest_point), flip_side(side)});
>>>>>>> dwn
      }
    }

    //TODO: now that we return multiple results we need to score
    //these lower or maybe not return them
    //if we have nothing because of heading we'll just ignore it
    if(correlated.size() == 0 && heading_filtered.size())
      for(auto& path_edge : heading_filtered)
        correlated.push_back(std::move(path_edge));
  }

  //give it back
  return correlated;
}

std::tuple<PointLL, float, size_t> project(const PointLL& p, const std::vector<PointLL>& shape) {
  size_t closest_segment = 0;
  float closest_distance = std::numeric_limits<float>::max();
  PointLL closest_point{};

  //for each segment
  for(size_t i = 0; i < shape.size() - 1; ++i) {
    //project a onto b where b is the origin vector representing this segment
    //and a is the origin vector to the point we are projecting, (a.b/b.b)*b
    const auto& u = shape[i];
    const auto& v = shape[i + 1];
    auto bx = v.first - u.first;
    auto by = v.second - u.second;
    auto sq = bx*bx + by*by;
    //avoid divided-by-zero which gives a NaN scale, otherwise comparisons below will fail
    const auto scale = sq > 0? (((p.first - u.first)*bx + (p.second - u.second)*by) / sq) : 0.f;
    //projects along the ray before u
    if(scale <= 0.f) {
      bx = u.first;
      by = u.second;
    }//projects along the ray after v
    else if(scale >= 1.f) {
      bx = v.first;
      by = v.second;
    }//projects along the ray between u and v
    else {
      bx = bx*scale + u.first;
      by = by*scale + u.second;
    }
    //check if this point is better
    PointLL point(bx, by);
    const auto distance = p.Distance(point);
    if(distance < closest_distance) {
      closest_segment = i;
      closest_distance = distance;
      closest_point = std::move(point);
    }
  }

  return std::make_tuple(std::move(closest_point), closest_distance, closest_segment);
}

//TODO: this is frought with peril. to properly to this we need to know
//where in the world we are and use lower casing rules that are appropriate
//so we can maximize the similarity measure.
//are the names similar enough to consider them matching
/*bool name_filter() {

}*/

// Test if this location is an isolated "island" without connectivity to the
// larger routing graph. Does a breadth first search - if possible paths are
// exhausted within some threshold this returns a set of edges within the
// island.
std::unordered_set<GraphId> island(const PathLocation& location,
             GraphReader& reader, const NodeFilter& node_filter,
             const EdgeFilter& edge_filter, const uint32_t edge_threshold,
             const uint32_t length_threshold, const uint32_t node_threshold) {
  std::unordered_set<GraphId> todo(edge_threshold);
  std::unordered_set<GraphId> done(edge_threshold);

  // Seed the list of edges to expand
  for (const auto& edge : location.edges) {
    todo.insert(edge.id);
  }

  // We are done if we hit a threshold meaning it isn't an island or we ran
  // out of edges and we determine it is an island
  uint32_t total_edge_length = 0;
  uint32_t nodes_expanded = 0;
  while ((done.size() < edge_threshold || total_edge_length < length_threshold ||
          nodes_expanded < node_threshold) && todo.size()) {
    // Get the next edge
    const GraphId edge = *todo.cbegin();
    done.emplace(edge);

    // Get the directed edge - filter it out if not accessible
    const DirectedEdge* directededge = reader.GetGraphTile(edge)->directededge(edge);
    if (edge_filter(directededge) == 0.0f) {
      continue;
    }
    total_edge_length += directededge->length();

    // Get the end node - filter it out if not accessible
    const GraphId node = directededge->endnode();
    const GraphTile* tile = reader.GetGraphTile(node);
    const NodeInfo* nodeinfo = tile->node(node);
    if (node_filter(nodeinfo)) {
      continue;
    }

    // Expand edges from the node
    bool expanded = false;
    GraphId edgeid(node.tileid(), node.level(), nodeinfo->edge_index());
    directededge = tile->directededge(nodeinfo->edge_index());
    for (uint32_t i = 0; i < nodeinfo->edge_count(); i++, directededge++, edgeid++) {
      // Skip transition edges, transit connection edges, and edges that are not allowed
      if (directededge->trans_up() || directededge->trans_down() ||
          directededge->use() == Use::kTransitConnection ||
          edge_filter(directededge) == 0.0f) {
        continue;
      }

      // Add to the todo list
      todo.emplace(edgeid);
      expanded = true;
    }
    nodes_expanded += expanded;
  }

  // If there are still edges to do then we broke out of the loop above due to
  // meeting thresholds and this is not a disconnected island. If there are no
  // more edges then this is a disconnected island and we want to know what
  // edges constitute the island so a second pass can avoid them
  return (todo.size() == 0) ? done : std::unordered_set<GraphId>{};
}

PathLocation search(const Location& location, GraphReader& reader,
    const EdgeFilter& edge_filter, const NodeFilter& node_filter) {
  //iterate over bins in a closest first manner
  const auto& tiles = reader.GetTileHierarchy().levels().rbegin()->second.tiles;
  const auto level = reader.GetTileHierarchy().levels().rbegin()->first;
  auto binner = tiles.ClosestFirst(location.latlng_);

  //TODO: if a name is supplied try to find edges with similar name
  //TODO: change filtering from boolean to floating point 0-1

  //give up if we find nothing after a while
  std::list<candidate_t> candidates;
  float best_so_far = std::numeric_limits<float>::max();
  //TODO: it might be more expensive to track this than not
  //std::unordered_set<uint64_t> seen(max_results * CONNECTED_REGION);
  size_t bins = 0;
  while(candidates.size() < MAX_CANDIDATES) {
    try {
      //TODO: make configurable the radius at which we give up searching
      auto bin = binner();
      if(std::get<2>(bin) > SEARCH_CUTOFF)
        break;

      //the closest thing in this bin is further than what we have already
      //TODO: keep searching a if we have very few or very crappy results so far
      if(candidates.size() && std::get<2>(bin) > best_so_far)
        break;

      //grab the tile the lat, lon is in
      const auto* tile = reader.GetGraphTile(GraphId(std::get<0>(bin), level, 0));
      if(!tile)
        continue;

      //iterate over the edges in the bin
      auto edges = tile->GetBin(std::get<1>(bin));
      bins += static_cast<size_t>(edges.size() > 0);
      for(auto e : edges) {
        //nothing new here
        //if(seen.find(e) != seen.cend())
        //  continue;
        //get the tile and edge
        candidate_t candidate;
        candidate.tile = e.tileid() != tile->id().tileid() ? reader.GetGraphTile(e) : tile;
        candidate.edge = candidate.tile->directededge(e);
        //no thanks on this one
        //seen.insert(e);
        if(edge_filter(candidate.edge) == 0.f) {
          //or its evil twin
          e = get_opposing(e, candidate.tile, reader);
          //seen.insert(e);
          if(!e.Is_Valid() || edge_filter(candidate.edge = candidate.tile->directededge(e)) == 0.f)
            continue;
        }
        //get some info about the edge
        candidate.id = e;
        candidate.edge_info = candidate.tile->edgeinfo(candidate.edge->edgeinfo_offset());
        candidate.point = project(location.latlng_, candidate.edge_info->shape());
        //add it in
        candidates.push_back(std::move(candidate));
        if(best_so_far > std::get<1>(candidate.point))
          best_so_far = std::get<1>(candidate.point);
      }
    }
    catch(...) {
      throw std::runtime_error("No data found for location");
    }
  }

  //keep track of bins we looked in but only the ones that had something
  //would rather log this in the service only, so lets figure a way to pass it back
  //midgard::logging::Log("valhalla_loki_bins_searched::" + std::to_string(bins), " [ANALYTICS] ");

  //keep the results and keep track of what edges are on islands
  std::list<std::unordered_set<uint64_t> > islands;
  auto on_island = [&islands](const GraphId& id) {
    for(const auto& i : islands)
      if(i.find(id) != i.cend())
        return true;
    return false;
  };

  //recursive depth first search for expanding edges
  std::function<void (const GraphTile*&, const DirectedEdge*, std::unordered_set<uint64_t>&)> depth_first;
  depth_first = [&depth_first, &reader, &node_filter](const GraphTile*& tile, const DirectedEdge* edge, std::unordered_set<uint64_t>& visited) {
      //can we get through the end node
      if(tile->id().tileid() != edge->endnode().tileid())
        tile = reader.GetGraphTile(edge->endnode());
      const auto* node = tile->node(edge->endnode());
      if(!node_filter(node)) {
        //for each edge
        iterable_t<const DirectedEdge> edges(tile->directededge(node->edge_index()), node->edge_count());
        for(const auto& e : edges){
          //saw this one yet?
          auto id = tile->id();
          id.fields.id = &e - tile->directededge(0);
          auto i = visited.insert(id);
          if(i.second && visited.size() < CONNECTED_REGION) {
            depth_first(tile, &e, visited);
          }
        }
      }
    };

  //check each candidate best first
  candidates.sort([](const candidate_t& a, const candidate_t& b){return std::get<1>(a.point) < std::get<1>(b.point);});
  PathLocation result(location);
  size_t result_count = 0;
  for(const auto& candidate : candidates) {
    //we've found enough results
    if(result_count - islands.size() > 0)
      break;

    //if this candidate was already found to be an island we dont really need it
    //since we have a better candidate from the same island earlier
    if(on_island(candidate.id))
      continue;

    //this may be at a node, either because it was the best thing or from snap tolerance
    bool front = std::get<0>(candidate.point) == candidate.edge_info->shape().front() ||
                 location.latlng_.Distance(candidate.edge_info->shape().front()) < NODE_SNAP;
    bool back = std::get<0>(candidate.point) == candidate.edge_info->shape().back() ||
                location.latlng_.Distance(candidate.edge_info->shape().back()) < NODE_SNAP;
    std::vector<PathLocation::PathEdge> correlated;
    //it was the begin node, so switch to opposing
    if((front && candidate.edge->forward()) || (back && !candidate.edge->forward())) {
      candidate_t opp_node(candidate, reader);
      opp_node.point = std::make_tuple(front ? candidate.edge_info->shape().front() : candidate.edge_info->shape().back(), 0, 0);
      std::get<1>(opp_node.point) = location.latlng_.Distance(std::get<0>(opp_node.point));
      correlated = correlate_node(reader, location, edge_filter, opp_node);
    }//it was the end node
    else if((back && candidate.edge->forward()) || (front && !candidate.edge->forward())) {
      candidate_t node(candidate);
      node.point = std::make_tuple(front ? candidate.edge_info->shape().front() : candidate.edge_info->shape().back(), 0, 0);
      std::get<1>(node.point) = location.latlng_.Distance(std::get<0>(node.point));
      correlated = correlate_node(reader, location, edge_filter, node);
    }//it was along the edge
    else {
      correlated = correlate_edge(reader, location, edge_filter, candidate);
    }

    //what did we get?
    if(!correlated.size())
      continue;
    result.edges.insert(result.edges.end(),
      std::make_move_iterator(correlated.begin()),
      std::make_move_iterator(correlated.end()));

    //do a depth first search keeping edges visited in case its an island
    std::unordered_set<uint64_t> island(CONNECTED_REGION);
    const auto* tile = candidate.tile;
    const auto* edge = candidate.edge;
    depth_first(tile, edge, island);
    //its a proper island
    if(island.size() < CONNECTED_REGION) {
      islands.emplace_back(std::move(island));
      result_count += 1;
    }//its not but lets mark the opposing so we dont use it again
    else {
      tile = candidate.tile;
      islands.emplace_back(std::unordered_set<uint64_t>{get_opposing(candidate.id, tile, reader)});
      result_count += 2;
    }
  }

  //if we still found nothing that is no good..
  if(result_count == 0)
    throw std::runtime_error("No suitable edges near location");

  return result;
}

}

namespace valhalla {
namespace loki {

baldr::PathLocation Search(const Location& location, GraphReader& reader,
    const EdgeFilter& edge_filter, const NodeFilter& node_filter) {
  //TODO: check if its an island and then search again
  return search(location, reader, edge_filter, node_filter);
}

}
}
