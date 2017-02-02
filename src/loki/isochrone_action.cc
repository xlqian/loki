#include "loki/service.h"
#include "loki/search.h"
#include <valhalla/baldr/datetime.h>
#include <boost/property_tree/json_parser.hpp>
#include <rapidjson/pointer.h>

using namespace prime_server;
using namespace valhalla::baldr;

namespace {
const headers_t::value_type CORS{"Access-Control-Allow-Origin", "*"};
}

namespace valhalla {
  namespace loki {

    void loki_worker_t::init_isochrones(const rapidjson::Document& request) {
      //strip off unused information
      parse_locations(request);
      if(locations.size() < 1)
        throw valhalla_exception_t{400, 120};
      for(auto& l : locations)
        l.heading_.reset();

      //make sure the isoline definitions are valid
      auto* contours = rapidjson::Pointer("/contours").Get(request);
      if(!contours)
        throw valhalla_exception_t{400, 113};
      //check that the number of contours is ok
      if(contours->GetArray().Size() > max_contours)
        throw valhalla_exception_t{400, 152, std::to_string(max_contours)};
      size_t prev = 0;
      for(const auto& contour : contours->GetArray()) {
        auto* time_ptr = rapidjson::Pointer("/time").Get(request);
        size_t c = time_ptr ? time_ptr->GetUint() : -1;
        if(c < prev || c == -1)
          throw valhalla_exception_t{400, 111};
        if(c > max_time)
          throw valhalla_exception_t{400, 151, std::to_string(max_time)};
        prev = c;
      }
      parse_costing(request);
    }

    worker_t::result_t loki_worker_t::isochrones(rapidjson::Document& request, http_request_info_t& request_info) {
      init_isochrones(request);
      //check that location size does not exceed max
      if (locations.size() > max_locations.find("isochrone")->second)
        throw valhalla_exception_t{400, 150, std::to_string(max_locations.find("isochrone")->second)};

      auto costing = request["costing"].GetString();
      auto* date_type =  rapidjson::Pointer("/date_time/type").Get(request);

      auto& allocator = request.GetAllocator();
      //default to current date_time for mm or transit.
      if (! date_type && (costing == "multimodal" || costing == "transit")) {
        date_type->Set(0);
      }

      //check the date stuff
      auto* date_time_value =rapidjson::Pointer("/date_time/value").Get(request);
      if (date_type) {
        //not yet on this
        if(date_type->GetInt() == 2) {
          jsonify_error({501, 142}, request_info);
        }
        //what kind
        switch(date_type->GetInt()) {
        case 0: //current
          rapidjson::GetValueByPointer(request, "/locations/0")->AddMember("date_time", "current", allocator);
          break;
        case 1: //depart
          if(!date_time_value)
            throw valhalla_exception_t{400, 160};
          if (!DateTime::is_iso_local(date_time_value->GetString()))
            throw valhalla_exception_t{400, 162};
          rapidjson::GetValueByPointer(request, "/locations/0")->AddMember("date_time", rapidjson::Value{*date_time_value, allocator}, allocator);
          break;
        default:
          throw valhalla_exception_t{400, 163};
          break;
        }
      }

      try{
        //correlate the various locations to the underlying graph
        const auto projections = loki::Search(locations, reader, edge_filter, node_filter);
        for(size_t i = 0; i < locations.size(); ++i) {
          auto path_tmp = ("/correlated_" + std::to_string(i)).c_str();
          rapidjson::Pointer(path_tmp).Set(request, projections.at(locations[i]).ToRapidJson(i, allocator));
        }
      }
      catch(const std::exception&) {
        throw valhalla_exception_t{400, 171};
      }

      //pass it on
      rapidjson::StringBuffer buffer;
      rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);
      request.Accept(writer);
      worker_t::result_t result{true};
      result.messages.emplace_back(buffer.GetString());

      return result;
    }

  }
}
