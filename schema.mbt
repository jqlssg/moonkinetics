///|
fn object_keys(
  value : Json,
  allowed : Array[String],
) -> Map[String, Json] raise KineticsError {
  let object = match value {
    Object(o) => o
    _ => raise Invalid("JSON object expected")
  }
  if object.length() != allowed.length() {
    raise Invalid("missing or extra JSON keys")
  }
  for key in object.keys() {
    if !allowed.contains(key) {
      raise Invalid("unknown JSON key")
    }
  }
  object
}

///|
fn required(
  object : Map[String, Json],
  key : String,
) -> Json raise KineticsError {
  match object.get(key) {
    Some(v) => v
    None => raise Invalid("required JSON key missing")
  }
}

///|
fn number(value : Json) -> Double raise KineticsError {
  match value {
    Number(x, ..) => {
      finite(x, "JSON number")
      x
    }
    _ => raise Invalid("JSON number expected")
  }
}

///|
fn array(value : Json) -> Array[Json] raise KineticsError {
  match value {
    Array(a) => a
    _ => raise Invalid("JSON array expected")
  }
}

///|
fn coefficients(value : Json) -> Array[Int32] raise KineticsError {
  let a = array(value)
  if a.is_empty() || a.length() > 64 {
    raise Invalid("coefficient vector size")
  }
  let result : Array[Int32] = []
  for entry in a {
    let x = number(entry)
    if x < 0.0 || x > 1000000.0 || x != x.to_int().to_double() {
      raise Invalid("coefficient must be a supported integer")
    }
    result.push(x.to_int())
  }
  result
}

///|
/// Versioned model schema: exactly schema/species/volume/reactions; each
/// channel exactly input/output/rate. Unknown keys and duplicate keys fail.
pub fn network_from_json(text : String) -> Network raise KineticsError {
  let o = object_keys(guarded_json(text), [
    "schema", "species", "volume", "reactions",
  ])
  if number(required(o, "schema")) != 1.0 {
    raise Invalid("unsupported model schema version")
  }
  let names : Array[String] = []
  for name in array(required(o, "species")) {
    match name {
      String(s) => names.push(s)
      _ => raise Invalid("species name must be a string")
    }
    if names.length() > 64 {
      raise Budget("JSON species budget")
    }
  }
  let reactions : Array[Reaction] = []
  for value in array(required(o, "reactions")) {
    let r = object_keys(value, ["input", "output", "rate"])
    reactions.push(
      Reaction::new(
        coefficients(required(r, "input")),
        coefficients(required(r, "output")),
        number(required(r, "rate")),
      ),
    )
    if reactions.length() > 1024 {
      raise Budget("JSON reaction budget")
    }
  }
  Network::new(names, reactions, volume=number(required(o, "volume")))
}

///|
pub fn Network::to_json(self : Network) -> Json {
  let channels = self.reactions.map(r => {
    let j : Json = { "input": r.input, "output": r.output, "rate": r.rate }
    j
  })
  {
    "schema": 1,
    "species": self.names,
    "volume": self.volume,
    "reactions": channels,
  }
}
