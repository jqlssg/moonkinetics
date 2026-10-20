///|
pub fn Trace::to_json(self : Trace) -> Json {
  { "schema": "trace-v1", "times": self.times, "values": self.values }
}

///|
/// Numeric long-form CSV avoids user-controlled spreadsheet formulas.
pub fn Trace::to_csv(self : Trace) -> String {
  let out = StringBuilder()
  out.write_string("time,species_index,value\n")
  for k in 0..<self.times.length() {
    for i in 0..<self.values[k].length() {
      out.write_string("\{self.times[k]},\{i},\{self.values[k][i]}\n")
    }
  }
  out.to_string()
}

///|
/// The event path stores initial counts plus channel/time events, not an
/// ambiguous interpolated table. Import revalidates and replays every event.
pub fn Path::to_json(self : Path) -> Json {
  {
    "schema": "path-v1",
    "initial": self.states[0],
    "times": self.times,
    "channels": self.channels,
    "horizon": self.horizon,
  }
}

///|
fn json_integer(value : Json, maximum : Int32) -> Int32 raise KineticsError {
  let x = number(value)
  if x < 0.0 || x > maximum.to_double() || x != x.to_int().to_double() {
    raise Invalid("JSON integer outside policy")
  }
  x.to_int()
}

///|
pub fn path_from_json(net : Network, text : String) -> Path raise KineticsError {
  let o = object_keys(guarded_json(text), [
    "schema", "initial", "times", "channels", "horizon",
  ])
  match required(o, "schema") {
    String("path-v1") => ()
    _ => raise Invalid("unsupported path schema")
  }
  let initial_json = array(required(o, "initial"))
  if initial_json.length() != net.names.length() {
    raise Invalid("path initial shape")
  }
  let initial : Array[Int32] = []
  for value in initial_json {
    initial.push(json_integer(value, 1000000000))
  }
  let time_json = array(required(o, "times"))
  let channel_json = array(required(o, "channels"))
  if time_json.length() > 100000 || time_json.length() != channel_json.length() {
    raise Invalid("path JSON event size")
  }
  let times : Array[Double] = []
  let channels : Array[Int32] = []
  for i in 0..<time_json.length() {
    times.push(number(time_json[i]))
    channels.push(json_integer(channel_json[i], 1023))
  }
  Path::from_events(
    net,
    initial,
    times,
    channels,
    number(required(o, "horizon")),
  )
}

///|
pub fn Path::to_csv(self : Path) -> String {
  let out = StringBuilder()
  out.write_string("event,time,channel,species_index,count\n")
  for k in 0..<self.states.length() {
    let time = if k == 0 { 0.0 } else { self.times[k - 1] }
    let channel = if k == 0 { -1 } else { self.channels[k - 1] }
    for i in 0..<self.states[k].length() {
      out.write_string("\{k},\{time},\{channel},\{i},\{self.states[k][i]}\n")
    }
  }
  out.to_string()
}

///|
pub fn Distribution::to_json(self : Distribution) -> Json {
  {
    "schema": "distribution-v1",
    "probabilities": self.probabilities,
    "poisson_tail_bound": self.tail,
    "terms": self.terms,
    "segments": self.segments,
  }
}
