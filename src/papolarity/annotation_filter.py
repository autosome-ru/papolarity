def parse_condition(condition_str):
    k,vs = condition_str.split('=', maxsplit=1)
    return (k, vs.split(','))

# Known issue: filters now treat all attribute values as strings
def create_record_filter(condition_config):
    k,vs = condition_config
    return lambda rec: str(rec.attributes.get(k, '')) in vs
