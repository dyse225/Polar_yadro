%% Симуляция полярного кода: SC, SCL и CA-SCL с разными размерами списков (PARFOR)
clear all;
close all;
clc;

% === ПАРАМЕТРЫ СИМУЛЯЦИИ ===
n_vec = 10;                           % Порядок блока (N = 2^n)
block_length = 2^n_vec;
info_length = 2^(n_vec - 1);
epsilon = 0.28;                       % Bhattacharyya параметр
ebno_vec = 0.5:0.25:5;              % Eb/N0 в дБ
crc_bits = 11;                        % CRC для CA-SCL

% Размеры списков для тестирования
list_sizes = [4, 8, 16];
num_lists = length(list_sizes);

% === ИНИЦИАЛИЗАЦИЯ ===
fprintf('Параметры симуляции:\n');
fprintf('  Размер блока N = %d\n', block_length);
fprintf('  Длина инфо K = %d\n', info_length);
fprintf('  Eb/N0 диапазон: [%.2f, %.2f] дБ\n\n', min(ebno_vec), max(ebno_vec));

% Запуск пула потоков (если еще не запущен)
if isempty(gcp('nocreate'))
    parpool; 
end

% Создание полярного кода (SC/SCL)
fprintf('Создание полярного кода...\n');
polar_code = Polar(block_length, info_length, epsilon, 0); 

% === ХРАНИЛИЩЕ РЕЗУЛЬТАТОВ ===
results = struct();
results.ebno = ebno_vec;

% === SC ДЕКОДЕР (L=1, без CRC) ===
% SC очень быстрый, parfor тут особо не нужен, оставим как есть или добавим, если ebno много
fprintf('SC декодер...\n');
[bler_sc, ber_sc] = polar_code.get_bler_quick(ebno_vec, 1);
results.sc.bler = bler_sc;
results.sc.ber = ber_sc;

% === SCL ДЕКОДЕРЫ (PARFOR) ===
fprintf('Запуск SCL декодеров в параллельном режиме...\n');

% Временные хранилища (cell array), так как parfor не любит прямую запись в struct(L)
scl_bler_temp = cell(1, num_lists);
scl_ber_temp  = cell(1, num_lists);

parfor idx = 1:num_lists
    L = list_sizes(idx);
    fprintf('  -> Воркер считает SCL L=%d...\n', L);
    % Важно: polar_code передается воркерам автоматически
    [b, e] = polar_code.get_bler_quick(ebno_vec, L);
    
    scl_bler_temp{idx} = b;
    scl_ber_temp{idx} = e;
end

% Перекладываем из временных ячеек в структуру results
for idx = 1:num_lists
    L = list_sizes(idx);
    results.scl(L).bler = scl_bler_temp{idx};
    results.scl(L).ber  = scl_ber_temp{idx};
end
fprintf('SCL готов.\n\n');

% === CA-SCL ДЕКОДЕРЫ (PARFOR) ===
fprintf('Создание CA-SCL кода с CRC=%d бит...\n', crc_bits);
polar_code_crc = PolarCode(block_length, info_length - crc_bits, epsilon, crc_bits);

fprintf('Запуск CA-SCL декодеров в параллельном режиме...\n');

cascl_bler_temp = cell(1, num_lists);
cascl_ber_temp  = cell(1, num_lists);

parfor idx = 1:num_lists
    L = list_sizes(idx);
    fprintf('  -> Воркер считает CA-SCL L=%d...\n', L);
    
    [b, e] = polar_code_crc.get_bler_quick(ebno_vec, L);
    
    cascl_bler_temp{idx} = b;
    cascl_ber_temp{idx} = e;
end

% Перекладываем результаты
for idx = 1:num_lists
    L = list_sizes(idx);
    results.cascl(L).bler = cascl_bler_temp{idx};
    results.cascl(L).ber  = cascl_ber_temp{idx};
end
fprintf('CA-SCL готов.\n');

% === ОЧИСТКА НУЛЕВЫХ ЗНАЧЕНИЙ ===
results.sc.ber(results.sc.ber == 0) = NaN;
for L = list_sizes
    results.scl(L).ber(results.scl(L).ber == 0) = NaN;
    results.cascl(L).ber(results.cascl(L).ber == 0) = NaN;
end

% === ПОСТРОЕНИЕ ГРАФИКА ===
fprintf('Построение графиков...\n');
p = figure('Position', [100, 100, 1000, 700]);
hold on;

% Цветовые схемы
colors_scl = [0, 0.8, 1; 0, 0.6, 0.8; 0, 0.4, 0.6; 0, 0.2, 0.4];
colors_cascl_bler = [1, 0, 0; 0.8, 0, 0; 0.6, 0, 0; 0.4, 0, 0];

% SC
plot(ebno_vec, results.sc.ber, 'Color', [0, 0.5, 1], 'LineStyle', '-', ...
    'Marker', 'o', 'MarkerSize', 6, 'LineWidth', 2, ...
    'DisplayName', 'SC');

% SCL
for idx = 1:length(list_sizes)
    L = list_sizes(idx);
    plot(ebno_vec, results.scl(L).ber, 'Color', colors_scl(idx, :), ...
        'LineStyle', '--', 'Marker', 's', 'MarkerSize', 5, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('SCL, L = %d', L));
end

% CA-SCL
for idx = 1:length(list_sizes)
    L = list_sizes(idx);
    plot(ebno_vec, results.cascl(L).ber, 'Color', colors_cascl_bler(idx, :), ...
        'LineStyle', '-', 'Marker', 'o', 'MarkerSize', 5, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('CA-SCL, L = %d', L));
end

hold off;

% Оформление
set(gca, 'YScale', 'log');
set(gca, 'FontSize', 11, 'FontName', 'Arial');
xlabel('E_b/N_0 (дБ)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('BER (Bit Error Rate)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Полярный код: N=%d, K=%d (с CRC=%d для CA-SCL)', ...
    block_length, info_length, crc_bits), 'FontSize', 13, 'FontWeight', 'bold');

leg = legend('Location', 'southwest', 'FontSize', 10);
leg.Box = 'on';
grid on;
set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.4);
xlim([min(ebno_vec) - 0.2, max(ebno_vec) + 0.2]);
ylim([1e-6, 1]);

saveas(p,'BER_SC_SCL_CA-SCL_Parfor.fig');
fprintf('Готово!\n');