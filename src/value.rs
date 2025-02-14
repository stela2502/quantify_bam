
/// Represents a single value in a bigWig file
#[derive(Copy, Clone, Debug, PartialEq)]
#[cfg_attr(feature = "write", derive(Serialize, Deserialize))]
pub struct Value {
    pub start: u32,
    pub end: u32,
    pub value: f32,
}

impl Value{
    pub fn flat(&self) -> ( u32, u32, f32) {
        ( self.start, self.end, self.value )
    }
}

