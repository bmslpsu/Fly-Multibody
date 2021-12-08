
#include <Servo.h>

#define daq_pin A1
#define servo_pin 9

int daq_signal;
float deg_val;

const float max_val = 92.2904; // set based on signal design
// const float max_val = 40; // set based on signal design

Servo myservo; // create servo object

void setup() {
  // 
  Serial.begin(9600);
  Serial.println("Initialize");
  myservo.attach(servo_pin); // attach servo
  myservo.write(0);
  delay(0.5);
}

void loop() {
  daq_signal = analogRead(daq_pin);
  deg_val = map(daq_signal, 0, 1023, 0, max_val);
  myservo.write(deg_val);
  Serial.println(deg_val);
  delay(0.1);

  Serial.print(deg_val);
  Serial.print("\n");
}
