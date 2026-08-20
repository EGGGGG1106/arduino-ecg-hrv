// ecg.ino - AD8232 ECG acquisition on Arduino UNO R4
//
// Wiring:
//   AD8232 OUTPUT -> A0
//   AD8232 LO+    -> D10
//   AD8232 LO-    -> D11
//
// Serial output:
//   banner : #ECG fs=500 bits=14
//   sample : <adc>,<leadOff>   e.g. "8192,0"

const uint32_t SAMPLE_RATE_HZ    = 500;
const uint32_t SAMPLE_INTERVAL_US = 1000000UL / SAMPLE_RATE_HZ;  // 2000 us
const uint8_t  ADC_BITS          = 14;

const uint8_t PIN_ECG      = A0;
const uint8_t PIN_LO_PLUS  = 10;
const uint8_t PIN_LO_MINUS = 11;

uint32_t nextSampleUs;

void setup() {
  Serial.begin(115200);

  // Native USB CDC: give the host a moment to open the port so the banner
  // is not lost, but do not block forever if nothing is listening.
  const uint32_t waitStart = millis();
  while (!Serial && (millis() - waitStart) < 3000) {
    ;
  }

  pinMode(PIN_LO_PLUS, INPUT);
  pinMode(PIN_LO_MINUS, INPUT);
  analogReadResolution(ADC_BITS);

  Serial.println("#ECG fs=500 bits=14");

  nextSampleUs = micros() + SAMPLE_INTERVAL_US;
}

void loop() {
  const uint32_t now = micros();

  // Signed difference keeps the comparison correct across micros() rollover.
  if ((int32_t)(now - nextSampleUs) < 0) {
    return;
  }

  nextSampleUs += SAMPLE_INTERVAL_US;

  // If we fell far behind (e.g. host stalled the USB write), resync instead
  // of bursting out a backlog of samples at the wrong timestamps.
  if ((int32_t)(now - nextSampleUs) > (int32_t)SAMPLE_INTERVAL_US) {
    nextSampleUs = now + SAMPLE_INTERVAL_US;
  }

  const uint16_t adc = analogRead(PIN_ECG);
  const uint8_t leadOff =
      (digitalRead(PIN_LO_PLUS) == HIGH || digitalRead(PIN_LO_MINUS) == HIGH) ? 1 : 0;

  Serial.print(adc);
  Serial.print(',');
  Serial.println(leadOff);
}
