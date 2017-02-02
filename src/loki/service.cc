#include <functional>
#include <string>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <cstdint>
#include <sstream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <valhalla/midgard/logging.h>
#include <valhalla/sif/autocost.h>
#include <valhalla/sif/bicyclecost.h>
#include <valhalla/sif/pedestriancost.h>
#include <valhalla/baldr/json.h>
#include <valhalla/baldr/errorcode_util.h>
#include <rapidjson/document.h>
#include <rapidjson/writer.h>
#include <rapidjson/pointer.h>

#include "loki/service.h"
#include "loki/search.h"

using namespace prime_server;
using namespace valhalla;
using namespace valhalla::midgard;
using namespace valhalla::baldr;
using namespace valhalla::sif;
using namespace valhalla::loki;

namespace {

  const std::unordered_map<std::string, loki_worker_t::ACTION_TYPE> PATH_TO_ACTION{
    {"/route", loki_worker_t::ROUTE},
    {"/viaroute", loki_worker_t::VIAROUTE},
    {"/locate", loki_worker_t::LOCATE},
    {"/one_to_many", loki_worker_t::ONE_TO_MANY},
    {"/many_to_one", loki_worker_t::MANY_TO_ONE},
    {"/many_to_many", loki_worker_t::MANY_TO_MANY},
    {"/sources_to_targets", loki_worker_t::SOURCES_TO_TARGETS},
    {"/optimized_route", loki_worker_t::OPTIMIZED_ROUTE},
    {"/isochrone", loki_worker_t::ISOCHRONE},
    {"/trace_route", loki_worker_t::TRACE_ROUTE},
    {"/trace_attributes", loki_worker_t::TRACE_ATTRIBUTES}
  };

  const headers_t::value_type CORS{"Access-Control-Allow-Origin", "*"};
  const headers_t::value_type JSON_MIME{"Content-type", "application/json;charset=utf-8"};
  const headers_t::value_type JS_MIME{"Content-type", "application/javascript;charset=utf-8"};

  rapidjson::Document from_request(const loki_worker_t::ACTION_TYPE& action, const http_request_t& request) {
    rapidjson::Document d;

    auto& allocator = d.GetAllocator();
    //parse the input
    const auto& json = request.query.find("json");
    if (json != request.query.end() && json->second.size()
      && json->second.front().size()) {
      std::istringstream is(json->second.front());
      d.Parse(is.str().c_str());
    }//no json parameter, check the body
    else if(!request.body.empty()) {
      std::istringstream is(request.body);
      d.Parse(is.str().c_str());
    }
    if (d.HasParseError())
      throw valhalla_exception_t{400, 100};

    // In case where the query is empty
    if (! d.IsObject()){
      d.SetObject();
    }

    //throw the query params into the ptree
    for(const auto& kv : request.query) {
      //skip json or empty entries
      if(kv.first == "json" || kv.first.size() == 0 || kv.second.size() == 0)
        continue;

      //turn single value entries into single key value
      if(kv.second.size() == 1) {
        d.AddMember(rapidjson::StringRef(kv.second.front().c_str()), rapidjson::StringRef(kv.first.c_str()), allocator);
        continue;
      }

      //make an array of values for this key
      rapidjson::Value array{rapidjson::kArrayType};
      for(const auto& value : kv.second) {
        array.PushBack(rapidjson::StringRef(value.c_str()), allocator);
      }
      d.AddMember(rapidjson::StringRef(kv.first.c_str()), array, allocator);
    }

    //if its osrm compatible lets make the location object conform to our standard input
    if(action == loki_worker_t::VIAROUTE) {
      auto& array = rapidjson::Pointer("/locations").Set(d, rapidjson::Value{rapidjson::kArrayType});

      for(const auto& location : d["loc"].GetArray()) {
        Location l = Location::FromCsv(location.GetString());
        rapidjson::Value ele{rapidjson::kObjectType};
        ele.AddMember("lon", l.latlng_.first, allocator)
            .AddMember("lat", l.latlng_.second, allocator);
        array.PushBack(ele, allocator);
      }
      d.RemoveMember("loc");
    }

    return d;
  }
}

namespace valhalla {
  namespace loki {
    worker_t::result_t loki_worker_t::jsonify_error(const valhalla_exception_t& exception, http_request_info_t& request_info) const {

      //build up the json map
      auto json_error = json::map({});
      json_error->emplace("status", exception.status_code_body);
      json_error->emplace("status_code", static_cast<uint64_t>(exception.status_code));
      json_error->emplace("error", std::string(exception.error_code_message));
      json_error->emplace("error_code", static_cast<uint64_t>(exception.error_code));

      //serialize it
      std::stringstream ss;
      if(jsonp)
        ss << *jsonp << '(';
      ss << *json_error;
      if(jsonp)
        ss << ')';

      worker_t::result_t result{false};
      http_response_t response(exception.status_code, exception.status_code_body, ss.str(), headers_t{CORS, jsonp ? JS_MIME : JSON_MIME});
      response.from_info(request_info);
      result.messages.emplace_back(response.to_string());

      return result;
    }

    void loki_worker_t::parse_locations(const rapidjson::Document& request) {
      //we require locations
      if (! request.HasMember("locations")){
        throw valhalla_exception_t{400, 110};
      }
      for(const auto& location : request["locations"].GetArray()) {
        try{
          locations.push_back(baldr::Location::FromRapidJson(location));
        }
        catch (...) {
          throw valhalla_exception_t{400, 130};
        }
      }
      valhalla::midgard::logging::Log("location_count::" + std::to_string(request["locations"].Size()), " [ANALYTICS] ");
    }

    void loki_worker_t::parse_costing(const rapidjson::Document& request) {
      //using the costing we can determine what type of edge filtering to use
      auto* consting_ptr = rapidjson::Pointer("/costing").Get(request);
      if (! consting_ptr) {
        throw valhalla_exception_t{400, 124};
      }
      auto costing = consting_ptr->GetString();
      valhalla::midgard::logging::Log(std::string("costing_type::")+ costing, " [ANALYTICS] ");
      // TODO - have a way of specifying mode at the location
      if(costing == "multimodal")
        costing = "pedestrian";

      // Get the costing options if in the config or get the empty default.
      // Creates the cost in the cost factory
      std::string method_options = std::string("/costing_options/") + costing;
      auto* method_options_ptr = rapidjson::Pointer(method_options.c_str()).Get(request);
      try{
        cost_ptr_t c;
        if (method_options_ptr){
          c = factory.Create(costing, *method_options_ptr);
        }else {
          c = factory.Create(costing, rapidjson::Value{});
        }

        edge_filter = c->GetEdgeFilter();
        node_filter = c->GetNodeFilter();
      }
      catch(const std::runtime_error&) {
        throw valhalla_exception_t{400, 125, std::string("'") + costing + "'"};
      }
    }

    loki_worker_t::loki_worker_t(const boost::property_tree::ptree& config):
        config(config), reader(config.get_child("mjolnir")), connectivity_map(config.get_child("mjolnir")),
        long_request(config.get<float>("loki.logging.long_request")),
        max_contours(config.get<unsigned int>("service_limits.isochrone.max_contours")),
        max_time(config.get<unsigned int>("service_limits.isochrone.max_time")),
        max_shape(config.get<size_t>("service_limits.trace.max_shape")) {

      // Keep a string noting which actions we support, throw if one isnt supported
      for (const auto& kv : config.get_child("loki.actions")) {
        auto path = "/" + kv.second.get_value<std::string>();
        if(PATH_TO_ACTION.find(path) == PATH_TO_ACTION.cend())
          throw valhalla_exception_t{400, 105, path};
        action_str.append("'" + path + "' ");
        actions.insert(path);
      }
      // Make sure we have at least something to support!
      if(action_str.empty())
        throw valhalla_exception_t{400, 102};

      //Build max_locations and max_distance maps
      for (const auto& kv : config.get_child("service_limits")) {
        if (kv.first != "skadi" && kv.first != "trace")
          max_locations.emplace(kv.first, config.get<size_t>("service_limits." + kv.first + ".max_locations"));
        if (kv.first != "skadi" && kv.first != "isochrone")
          max_distance.emplace(kv.first, config.get<float>("service_limits." + kv.first + ".max_distance"));
      }
      //this should never happen
      if (max_locations.empty())
        throw valhalla_exception_t{400, 103};

      if (max_distance.empty())
        throw valhalla_exception_t{400, 104};

      min_transit_walking_dis =
        config.get<int>("service_limits.pedestrian.min_transit_walking_distance");
      max_transit_walking_dis =
        config.get<int>("service_limits.pedestrian.max_transit_walking_distance");

      // Register edge/node costing methods
      // TODO: move this into the loop above
      factory.Register("auto", sif::CreateAutoCost);
      factory.Register("auto_shorter", sif::CreateAutoShorterCost);
      factory.Register("bus", sif::CreateBusCost);
      factory.Register("bicycle", sif::CreateBicycleCost);
      factory.Register("hov", sif::CreateHOVCost);
      factory.Register("pedestrian", sif::CreatePedestrianCost);
      factory.Register("truck", sif::CreateTruckCost);
      factory.Register("transit", sif::CreateTransitCost);

    }

    worker_t::result_t loki_worker_t::work(const std::list<zmq::message_t>& job, void* request_info, const worker_t::interrupt_function_t&) {
      //get time for start of request
      auto s = std::chrono::system_clock::now();
      auto& info = *static_cast<http_request_info_t*>(request_info);
      LOG_INFO("Got Loki Request " + std::to_string(info.id));

      try{
        //request parsing
        auto request = http_request_t::from_string(static_cast<const char*>(job.front().data()), job.front().size());

        //block all but get and post
        if(request.method != method_t::POST && request.method != method_t::GET)
          return jsonify_error({405, 101}, info);

        //is the request path action in the action set?
        auto action = PATH_TO_ACTION.find(request.path);
        if (action == PATH_TO_ACTION.cend() || actions.find(request.path) == actions.cend())
          return jsonify_error({404, 106, action_str}, info);

        //parse the query's json
        auto request_rj = from_request(action->second, request);

        auto* jsonp_ptr = rapidjson::Pointer("/jsonp").Get(request_rj);
        jsonp = jsonp_ptr ? boost::optional<std::string>{jsonp_ptr->GetString()} : boost::optional<std::string>{};
        //let further processes more easily know what kind of request it was
        rapidjson::SetValueByPointer(request_rj, "/action", action->second);
        //let further processes know about tracking
        auto do_not_track = request.headers.find("DNT");
        info.spare = do_not_track != request.headers.cend() && do_not_track->second == "1";

        worker_t::result_t result{false};
        //do request specific processing
        switch (action->second) {
          case ROUTE:
          case VIAROUTE:
            result = route(request_rj, info);
            break;
          case LOCATE:
            result = locate(request_rj, info);
            break;
          case ONE_TO_MANY:
          case MANY_TO_ONE:
          case MANY_TO_MANY:
          case SOURCES_TO_TARGETS:
          case OPTIMIZED_ROUTE:
            result = matrix(action->second, request_rj, info);
            break;
          case ISOCHRONE:
            result = isochrones(request_rj, info);
            break;
          case TRACE_ATTRIBUTES:
          case TRACE_ROUTE:
            result = trace_route(request_rj, info);
            break;
          default:
            //apparently you wanted something that we figured we'd support but havent written yet
            return jsonify_error({501, 107}, info);
        }
        //get processing time for loki
        auto e = std::chrono::system_clock::now();
        std::chrono::duration<float, std::milli> elapsed_time = e - s;
        //log request if greater than X (ms)
        auto work_units = locations.size() ? locations.size() : 1;
        if (!info.spare && elapsed_time.count() / work_units > long_request) {
          rapidjson::StringBuffer buffer;
          rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);
          request_rj.Accept(writer);
          LOG_WARN("loki::request elapsed time (ms)::"+ std::to_string(elapsed_time.count()));
          LOG_WARN(std::string("loki::request exceeded threshold::")+ buffer.GetString());
          midgard::logging::Log("valhalla_loki_long_request", " [ANALYTICS] ");
        }

        return result;
      }
      catch(const valhalla_exception_t& e) {
        valhalla::midgard::logging::Log("400::" + std::string(e.what()), " [ANALYTICS] ");
        return jsonify_error({e.status_code, e.error_code, e.extra}, info);
      }
      catch(const std::exception& e) {
        valhalla::midgard::logging::Log("400::" + std::string(e.what()), " [ANALYTICS] ");
        return jsonify_error({400, 199, std::string(e.what())}, info);
      }
    }

    void loki_worker_t::cleanup() {
      jsonp = boost::none;
      locations.clear();
      sources.clear();
      targets.clear();
      shape.clear();
      if(reader.OverCommitted())
        reader.Clear();
    }

    void run_service(const boost::property_tree::ptree& config) {
      //gets requests from the http server
      auto upstream_endpoint = config.get<std::string>("loki.service.proxy") + "_out";
      //sends them on to thor
      auto downstream_endpoint = config.get<std::string>("thor.service.proxy") + "_in";
      //or returns just location information back to the server
      auto loopback_endpoint = config.get<std::string>("httpd.service.loopback");
      auto interrupt_endpoint = config.get<std::string>("httpd.service.interrupt");

      //listen for requests
      zmq::context_t context;
      loki_worker_t loki_worker(config);
      prime_server::worker_t worker(context, upstream_endpoint, downstream_endpoint, loopback_endpoint, interrupt_endpoint,
        std::bind(&loki_worker_t::work, std::ref(loki_worker), std::placeholders::_1, std::placeholders::_2, std::placeholders::_3),
        std::bind(&loki_worker_t::cleanup, std::ref(loki_worker)));
      worker.work();

      //TODO: should we listen for SIGINT and terminate gracefully/exit(0)?
    }
  }
}
