///|
/// Bounded prescan before the core recursive parser. Duplicate object keys
/// are compared after JSON escape decoding, including escaped spellings.
fn guarded_json(text : String) -> Json raise KineticsError {
  if text.length() > 1000000 {
    raise Budget("model JSON exceeds one million UTF-16 units")
  }
  let stack : Array[Map[String, Bool]?] = []
  let mut quoted = false
  let mut escaped = false
  let mut token = StringBuilder()
  let mut pending : String? = None
  for ch in text.iter() {
    if quoted {
      token.write_char(ch)
      if escaped {
        escaped = false
      } else if ch == '\\' {
        escaped = true
      } else if ch == '"' {
        quoted = false
        let decoded = @json.parse(token.to_string()) catch {
          _ => raise Invalid("invalid JSON string")
        }
        pending = match decoded {
          String(s) => Some(s)
          _ => None
        }
      }
      continue
    }
    if ch == ' ' || ch == '\n' || ch == '\r' || ch == '\t' {
      continue
    }
    if ch == ':' {
      if stack.is_empty() {
        raise Invalid("JSON colon without object")
      }
      match (stack[stack.length() - 1], pending) {
        (Some(keys), Some(key)) => {
          if keys.contains(key) {
            raise Invalid("duplicate JSON object key")
          }
          keys[key] = true
        }
        _ => raise Invalid("JSON object key expected")
      }
      pending = None
    } else if ch == '"' {
      token = StringBuilder()
      token.write_char(ch)
      quoted = true
      pending = None
    } else {
      pending = None
      if ch == '{' {
        stack.push(Some(Map([])))
      } else if ch == '[' {
        stack.push(None)
      } else if ch == '}' || ch == ']' {
        if stack.is_empty() {
          raise Invalid("unbalanced JSON")
        }
        let _ = stack.pop()
      }
      if stack.length() > 6 {
        raise Budget("model JSON nesting exceeds six")
      }
    }
  }
  @json.parse(text) catch {
    _ => raise Invalid("invalid model JSON")
  }
}
